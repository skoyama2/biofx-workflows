version 1.0

import "BatchCramFiles.wdl" as batch_files

workflow Glimpse2LowPassImputationCrams {

    input {

        String pipeline_version = "0.0.1"

        # TAB-separated: cram_path<TAB>crai_path<TAB>sampleID (one sample per line).
        # Optional header line; set manifest_has_header accordingly.
        File sample_manifest

        File reference_chunks
        String output_basename

        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict

        String docker

        Boolean manifest_has_header = true

        # When set > 0 and sample count exceeds this value, samples are split into
        # ceil(N / max_files_per_ligate) batches. Each batch runs phase per reference
        # chunk, then one ligate per batch (output_basename.batch0, .batch1, ...).
        # 0 or negative = no sample batching.
        Int max_files_per_ligate = 0

        # If true, merge per-batch VCFs into one with bcftools merge (after all batches).
        Boolean merge_batch_outputs = false

    }

    Int header_offset = if manifest_has_header then 1 else 0
    Int total_lines = length(read_lines(sample_manifest))
    Int total_crams = total_lines - header_offset
    Int batch_size = if max_files_per_ligate > 0 then max_files_per_ligate else total_crams
    Int num_batches = if max_files_per_ligate <= 0 || total_crams <= max_files_per_ligate then 1 else (total_crams + max_files_per_ligate - 1) / max_files_per_ligate

    Array[Array[String]] manifest_rows = read_tsv(sample_manifest)

    scatter (batch_idx in range(num_batches)) {

        call batch_files.BatchCramFiles {
            input:
                manifest_rows = manifest_rows,
                header_offset = header_offset,
                batch_idx = batch_idx,
                batch_size = batch_size,
                total_crams = total_crams
        }

        scatter (reference_chunk in read_lines(reference_chunks)) {

            call GlimpsePhase {
                input:
                    crams = BatchCramFiles.crams,
                    crais = BatchCramFiles.crais,
                    reference_chunk = reference_chunk,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_fasta_dict = ref_fasta_dict,
                    docker = docker
            }

        }

        call GlimpseLigate {
            input:
                imputed_chunks = GlimpsePhase.imputed_vcf,
                imputed_chunks_indices = GlimpsePhase.imputed_vcf_index,
                ref_fasta_dict = ref_fasta_dict,
                output_basename = output_basename + ".batch" + batch_idx,
                docker = docker
        }

    }

    if (merge_batch_outputs && num_batches > 1) {
        call MergeBatchVcfs {
            input:
                batch_vcfs = GlimpseLigate.imputed_vcf,
                output_basename = output_basename,
                ref_fasta_dict = ref_fasta_dict,
                docker = docker
        }
    }

    output {
        Array[File] batch_imputed_vcf = GlimpseLigate.imputed_vcf
        Array[File] batch_imputed_vcf_index = GlimpseLigate.imputed_vcf_index
        File? merged_imputed_vcf = MergeBatchVcfs.merged_vcf
        File? merged_imputed_vcf_index = MergeBatchVcfs.merged_vcf_index
    }

}

task MergeBatchVcfs {

    input {
        Array[File] batch_vcfs
        String output_basename
        File ref_fasta_dict
        String docker
        Int mem_gb = 16
        Int disk_size_gb = 200
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 2
    }

    command <<<

        set -euo pipefail
        bcftools merge -l ~{write_lines(batch_vcfs)} -Oz -o merged_raw.vcf.gz
        bcftools sort merged_raw.vcf.gz -Ou \
            | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou \
            | bcftools norm -d both -Oz -o merged_cleaned.vcf.gz
        bcftools view -h --no-version merged_cleaned.vcf.gz > old_header.vcf
        java -jar /picard.jar UpdateVcfSequenceDictionary \
            -I old_header.vcf \
            --SD ~{ref_fasta_dict} -O new_header.vcf
        bcftools reheader -h new_header.vcf -o ~{output_basename}.merged.imputed.vcf.gz merged_cleaned.vcf.gz
        tabix ~{output_basename}.merged.imputed.vcf.gz

    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        File merged_vcf = "~{output_basename}.merged.imputed.vcf.gz"
        File merged_vcf_index = "~{output_basename}.merged.imputed.vcf.gz.tbi"
    }

}

task GlimpseLigate {

    input {

        String output_basename

        File ref_fasta_dict

        Array[File] imputed_chunks
        Array[File] imputed_chunks_indices

        Int mem_gb = 16
        Int disk_size_gb = 100
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 2

        String glimpse_ligate = "/usr/local/bin/GLIMPSE2_ligate"
        String docker

    }

    command <<<

    ~{glimpse_ligate}  \
        --input ~{write_lines(imputed_chunks)} \
        --output ligated.vcf.gz

    bcftools sort ligated.vcf.gz -Ou \
        | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou \
        | bcftools norm -d both -Oz -o ligated_cleaned.vcf.gz

    bcftools view \
        -h --no-version ligated_cleaned.vcf.gz > old_header.vcf

    java -jar /picard.jar UpdateVcfSequenceDictionary \
        -I old_header.vcf  \
        --SD ~{ref_fasta_dict} -O new_header.vcf

    bcftools reheader \
        -h new_header.vcf \
        -o ~{output_basename}.imputed.vcf.gz \
        ligated_cleaned.vcf.gz

    tabix ~{output_basename}.imputed.vcf.gz

    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        File imputed_vcf = "~{output_basename}.imputed.vcf.gz"
        File imputed_vcf_index = "~{output_basename}.imputed.vcf.gz.tbi"
    }

}

task GlimpsePhase {

    input {

        File reference_chunk

        Array[File] crams
        Array[File] crais

        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict

        Int mem_gb = 16
        Int disk_size_gb = 100
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 2

        String docker
        String glimpse_phase = "/usr/local/bin/GLIMPSE2_phase"

    }

    command <<<

        set -euo pipefail
        sort -V ~{write_lines(crams)} > sorted_crams.list

        ~{glimpse_phase}  \
            --bam-list sorted_crams.list \
            --reference ~{reference_chunk} \
            --output phase_output.bcf \
            --threads 1

    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        File imputed_vcf = "phase_output.bcf"
        File imputed_vcf_index = "phase_output.bcf.csi"
    }

}
