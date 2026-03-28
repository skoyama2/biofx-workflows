version 1.0

workflow BuildGlimpseReferenceChunks {
  input {
    Array[String] chromosomes = [
      "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
      "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
      "chrX"
    ]

    Array[File] reference_bcfs
    Array[File] reference_csis
    Array[File] map_files

    Float window_mb = 5.0
    Float window_cm = 5.0
    Int window_count = 150000

    String chrX_nonpar_region = "chrX:2781480-155701382"
  }

  scatter (i in range(length(chromosomes))) {
    call BuildOneChromosomeReference {
      input:
        chr = chromosomes[i],
        reference_bcf = reference_bcfs[i],
        reference_csi = reference_csis[i],
        map_file = map_files[i],
        window_mb = window_mb,
        window_cm = window_cm,
        window_count = window_count,
        chunk_region = if chromosomes[i] == "chrX" then chrX_nonpar_region else chromosomes[i]
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
    String chunk_region

    File reference_bcf
    File reference_csi
    File map_file

    Float window_mb
    Float window_cm
    Int window_count

    String docker = "skoyamamd/glimpse2:latest"
    String glimpse_chunk = "/usr/local/bin/GLIMPSE2_chunk"
    String glimpse_split_reference = "/usr/local/bin/GLIMPSE2_split_reference"

    Int disk_gb = 100
    Int memory_gb = 16
    Int cpu = 4
  }

  String vcf_filename = basename(reference_bcf)
  String vcf_csi_filename = basename(reference_csi)
  String vcf_stem = sub(vcf_filename, "\\.bcf$", "")

  command <<<
    set -euo pipefail

    mkdir -p phased_reference binary_reference_panel

    ln -s "~{reference_bcf}" "phased_reference/~{vcf_filename}"
    ln -s "~{reference_csi}" "phased_reference/~{vcf_csi_filename}"

    ~{glimpse_chunk} \
      --input "phased_reference/~{vcf_filename}" \
      --map "~{map_file}" \
      --sequential \
      --region "~{chunk_region}" \
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
    File chunk_txt = "binary_reference_panel/" + vcf_stem + ".chunk.txt"
    Array[File] split_outputs = glob("binary_reference_panel/" + vcf_stem + "*")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
  }

}
