version 1.0

workflow BuildGlimpseReferenceChunks {

  input {

    Array[String] chromosomes = [
      "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
      "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
    ]

    # VCF path is constructed as: vcf_prefix + chr + vcf_suffix
    # e.g. "gs://bucket/hgdp1kgp_" + "chr1" + ".filtered.SNV_INDEL.phased.shapeit5.bcf"
    String vcf_prefix = "gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_"
    String vcf_suffix = ".filtered.SNV_INDEL.phased.shapeit5.bcf"

    # One genetic map file per chromosome, in the same order as chromosomes
    Array[File] map_files

    Float window_mb = 5.0
    Float window_cm = 5.0
    Int window_count = 150000

  }

  scatter (i in range(length(chromosomes))) {
    call BuildOneChromosomeReference {
      input:
        chr = chromosomes[i],
        vcf_prefix = vcf_prefix,
        vcf_suffix = vcf_suffix,
        map_file = map_files[i],
        window_mb = window_mb,
        window_cm = window_cm,
        window_count = window_count
    }
  }

  output {
    Array[File] chunk_txts = BuildOneChromosomeReference.chunk_txt
    Array[Array[File]] split_outputs = BuildOneChromosomeReference.split_outputs
  }

}

task BuildOneChromosomeReference {
  input {
    String chr

    String vcf_prefix
    String vcf_suffix

    File map_file

    Float window_mb
    Float window_cm
    Int window_count

    String docker_image = "us.gcr.io/broad-dsde-methods/glimpse:kachulis_ck_bam_reader_retry_cf5822c"

    String glimpse_chunk = "/bin/GLIMPSE_chunk"
    String glimpse_split_reference = "/bin/GLIMPSE_split_reference"

    Int disk_gb = 100
    Int memory_gb = 16
    Int cpu = 4

  }

  String vcf_gcs_path = vcf_prefix + chr + vcf_suffix
  String vcf_filename = basename(vcf_gcs_path)
  String vcf_stem = basename(vcf_gcs_path, ".bcf")

  command <<<
    set -euo pipefail

    mkdir -p phased_reference binary_reference_panel

    gcloud storage cp -n "~{vcf_gcs_path}" phased_reference/
    gcloud storage cp -n "~{vcf_gcs_path}.csi" phased_reference/

    ~{glimpse_chunk} \
      --input "phased_reference/~{vcf_filename}" \
      --map "~{map_file}" \
      --sequential \
      --region "~{chr}" \
      --window-mb ~{window_mb} \
      --window-cm ~{window_cm} \
      --window-count ~{window_count} \
      --output "binary_reference_panel/~{vcf_stem}.chunk.txt"

    cut -f 2,3,4 "binary_reference_panel/~{vcf_stem}.chunk.txt" \
      | while read input_chr buf imp; do

        ~{glimpse_split_reference} \
          --reference "phased_reference/~{vcf_filename}" \
          --map "~{map_file}" \
          --input-region "${buf}" \
          --output-region "${imp}" \
          --threads 1 \
          --output "binary_reference_panel/~{vcf_stem}"
      done
  >>>

  output {
    File chunk_txt = "binary_reference_panel/~{vcf_stem}.chunk.txt"
    Array[File] split_outputs = glob("binary_reference_panel/~{vcf_stem}*")
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
  }
}
