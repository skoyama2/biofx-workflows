version 1.0

import "BatchCramFiles.wdl" as batch_files

workflow Glimpse2LowPassImputationCrams {

    input {

        String pipeline_version = "0.0.1"

        # TAB-separated (one sample per line): col1 cram, col2 crai, col3 sample name (optional).
        # When col3 is present (non-empty), it is written as column 2 of GLIMPSE2_phase --bam-list.
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
        # chunk, then one ligate per batch (output_basename.batch0.imputed.bcf, ...).
        # 0 or negative = no sample batching: single ligate uses output_basename.imputed.bcf
        # and is exposed as merged_imputed_bcf (not batch_imputed_bcf).

        Int max_files_per_ligate = 0

        # If true, merge per-batch BCFs into one with bcftools merge (after all batches).

        Boolean merge_batch_outputs = false

        # Extra attempts after failure (preemptible, etc.). On Terra (Google backend), pair with
        # workflow options JSON memory_retry_multiplier (see Glimpse2LowPassImputationCrams.workflow_options.json)
        # so stderr OOM patterns trigger a retry with multiplied memory.

        Int phase_max_retries = 1
        Int ligate_max_retries = 2
        Int merge_batch_max_retries = 1

    }

    Int header_offset = if manifest_has_header then 1 else 0

    # Must match BatchCramFiles indexing (read_tsv rows). read_lines(sample_manifest) can disagree
    # (e.g. blank lines) and made only the last batch line up — batch0/1 then fail, last batch succeeds.

    Array[Array[String]] manifest_rows = read_tsv(sample_manifest)
    Int manifest_nrows = length(manifest_rows)
    Int total_crams = manifest_nrows - header_offset
    Int batch_size = if max_files_per_ligate > 0 then max_files_per_ligate else total_crams
    Int num_batches = if max_files_per_ligate <= 0 || total_crams <= max_files_per_ligate then 1 else (total_crams + max_files_per_ligate - 1) / max_files_per_ligate

    # Hoist once: avoids Terra/Cromwell oddities with read_lines nested inside the batch scatter.
    Array[String] reference_chunk_lines = read_lines(reference_chunks)

    scatter (batch_idx in range(num_batches)) {

        call batch_files.BatchCramFiles {
            input:
                manifest_rows = manifest_rows,
                header_offset = header_offset,
                batch_idx = batch_idx,
                batch_size = batch_size,
                total_crams = total_crams
        }

        scatter (reference_chunk in reference_chunk_lines) {

            call GlimpsePhase {
                input:
                    crams = BatchCramFiles.crams,
                    crais = BatchCramFiles.crais,
                    cram_phase_sample_names = BatchCramFiles.phase_sample_names,
                    reference_chunk = reference_chunk,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_fasta_dict = ref_fasta_dict,
                    max_retries = phase_max_retries,
                    docker = docker
            }

        }

        call GlimpseLigate {
            input:
                imputed_chunks = GlimpsePhase.imputed_bcf,
                imputed_chunks_indices = GlimpsePhase.imputed_bcf_index,
                ref_fasta_dict = ref_fasta_dict,
                output_basename = if num_batches > 1 then output_basename + ".batch" + batch_idx else output_basename,
                max_retries = ligate_max_retries,
                docker = docker
        }

    }

    if (merge_batch_outputs && num_batches > 1) {
        call MergeBatchBcfs {
            input:
                batch_bcfs = GlimpseLigate.imputed_bcf,
                batch_bcf_indices = GlimpseLigate.imputed_bcf_index,
                output_basename = output_basename,
                ref_fasta_dict = ref_fasta_dict,
                max_retries = merge_batch_max_retries,
                docker = docker
        }
    }

    output {
        Array[File]? batch_bcf = if num_batches > 1 then GlimpseLigate.imputed_bcf else []
        Array[File]? batch_bcf_index = if num_batches > 1 then GlimpseLigate.imputed_bcf_index else []
        File? imputed_bcf = if num_batches == 1 then GlimpseLigate.imputed_bcf[0] else MergeBatchBcfs.merged_bcf
        File? imputed_bcf_index = if num_batches == 1 then GlimpseLigate.imputed_bcf_index[0] else MergeBatchBcfs.merged_bcf_index
    }

}

task MergeBatchBcfs {

    input {
        Array[File] batch_bcfs
        # Same gather order as batch_bcfs (per-batch ligate outputs).
        Array[File] batch_bcf_indices
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
        bcf_list=~{write_lines(batch_bcfs)}
        csi_list=~{write_lines(batch_bcf_indices)}
        n_bcf=$(wc -l < "${bcf_list}" | tr -d ' ')
        n_csi=$(wc -l < "${csi_list}" | tr -d ' ')
        if [ "${n_bcf}" != "${n_csi}" ]; then
            echo "ERROR: batch BCF count (${n_bcf}) != index count (${n_csi})" >&2
            exit 1
        fi

        mkdir -p merge_staging
        merge_paths=merge.bcf.paths

        : > "${merge_paths}"

        # Pair paths without process substitution (sh-safe): paste to a file, then read in this shell.

        paste "${bcf_list}" "${csi_list}" > bcf_csi.pairs

        i=0

        while IFS="$(printf '\t')" read -r bcf csi || [ -n "${bcf}" ]; do
            bcf=$(printf '%s' "${bcf}" | tr -d '\r')
            csi=$(printf '%s' "${csi}" | tr -d '\r')
            [ -n "${bcf}" ] && [ -n "${csi}" ] || continue
            out="merge_staging/batch_${i}.bcf"
            cp -f "${bcf}" "${out}"
            cp -f "${csi}" "${out}.csi"
            echo "${PWD}/${out}" >> "${merge_paths}"
            i=$((i + 1))
        done < bcf_csi.pairs

        test -s "${merge_paths}" || { echo "ERROR: no batch BCFs to merge" >&2; exit 1; }

        bcftools merge -l "${merge_paths}" -Ob -o merged_raw.bcf

        bcftools sort merged_raw.bcf -Ou \
            | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou \
            | bcftools norm -d both -Ob -o merged_cleaned.bcf

        bcftools view -h --no-version merged_cleaned.bcf > old_header.vcf

        java -jar /picard.jar UpdateVcfSequenceDictionary \
            -I old_header.vcf \
            --SD ~{ref_fasta_dict} \
            -O new_header.vcf

        bcftools reheader -h new_header.vcf -o ~{output_basename}.imputed.bcf merged_cleaned.bcf

        bcftools index -f ~{output_basename}.imputed.bcf

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
        File merged_bcf = "~{output_basename}.imputed.bcf"
        File merged_bcf_index = "~{output_basename}.imputed.bcf.csi"
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

    set -euo pipefail
    trap 'echo "ERROR: GlimpseLigate failed at line ${LINENO}" >&2' ERR
    chunks_raw_list=~{write_lines(imputed_chunks)}
    chunks_sorted_list=chunks.sorted.list
    chunks_meta=chunks.meta.tsv
    echo "INFO: raw chunk count = $(wc -l < "${chunks_raw_list}")"

    # Keep ligate input order genomic (not call/shard order), otherwise GLIMPSE2 can fail with:
    # "Overlap is empty".

    : > "${chunks_meta}"

    while IFS= read -r bcf; do
        [ -n "${bcf}" ] || continue
        echo "INFO: scan chunk ${bcf}"
        first_record=$(bcftools query -f '%CHROM\t%POS\n' "${bcf}" | awk 'NR==1 {print; exit}' || true)
        chrom=$(printf '%s\n' "${first_record}" | cut -f1)
        pos=$(printf '%s\n' "${first_record}" | cut -f2)
        if [ -z "${chrom}" ] || [ -z "${pos}" ]; then
            echo "ERROR: empty or unreadable chunk BCF: ${bcf}" >&2
            exit 1
        fi
        printf '%s\t%s\t%s\n' "${chrom}" "${pos}" "${bcf}" >> "${chunks_meta}"
    done < "${chunks_raw_list}"

    sort -k1,1V -k2,2n "${chunks_meta}" | cut -f3 > "${chunks_sorted_list}"
    echo "INFO: sorted chunk list"
    cat "${chunks_sorted_list}"

    # GLIMPSE2_ligate writes BCF (same as phase inputs); do not use a .vcf.gz name or bcftools mis-detects format.

    ~{glimpse_ligate}  \
        --input "${chunks_sorted_list}" \
        --output ligated.bcf

    bcftools sort ligated.bcf -Ou \
        | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou \
        | bcftools norm -d both -Ob -o ligated_cleaned.bcf

    bcftools view \
        -h --no-version ligated_cleaned.bcf > old_header.vcf

    java -jar /picard.jar UpdateVcfSequenceDictionary \
        -I old_header.vcf  \
        --SD ~{ref_fasta_dict} -O new_header.vcf

    bcftools reheader \
        -h new_header.vcf \
        -o ~{output_basename}.imputed.bcf \
        ligated_cleaned.bcf

    bcftools index -f ~{output_basename}.imputed.bcf

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
        File imputed_bcf = "~{output_basename}.imputed.bcf"
        File imputed_bcf_index = "~{output_basename}.imputed.bcf.csi"
    }

}

task GlimpsePhase {

    input {

        File reference_chunk

        Array[File] crams
        Array[File] crais
        # Parallel to crams; from manifest col3 when present. Empty string => single-column bam-list line.
        Array[String] cram_phase_sample_names

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

        cram_list=~{write_lines(crams)}
        crai_list=~{write_lines(crais)}
        sample_list=~{write_lines(cram_phase_sample_names)}

        mkdir -p staged_crams
        
        while read -r cram_path; do
          ln -sf "${cram_path}" "staged_crams/$(basename "${cram_path}")"
        done < "${cram_list}"

        while read -r crai_path; do
          ln -sf "${crai_path}" "staged_crams/$(basename "${crai_path}")"
        done < "${crai_list}"

        # Pair cram paths with optional sample IDs, sort by staged CRAM basename (same as sort -V on paths).

        paste "${cram_list}" "${sample_list}" \
          | while IFS="$(printf '\t')" read -r cram_path sample || [ -n "${cram_path}" ]; do
          
          # If cram_path is empty, skip the line
          if [ -z "${cram_path}" ]; then
            continue
          fi

          staged="${PWD}/staged_crams/$(basename "${cram_path}")"
          base="$(basename "${cram_path}")"

          if [ -n "${sample}" ]; then
            printf '%s\t%s\t%s\n' "${base}" "${staged}" "${sample}"
          else
            printf '%s\t%s\n' "${base}" "${staged}"
          fi

        done > bam_list.unsorted.tsv

        sort -t "$(printf '\t')" -k1,1V bam_list.unsorted.tsv \
          | while IFS="$(printf '\t')" read -r _ staged sample_rest || [ -n "${staged}" ]; do

          # If staged is empty, skip the line
          if [ -z "${staged}" ]; then
            continue
          fi

          if [ -n "${sample_rest}" ]; then
            printf '%s\t%s\n' "${staged}" "${sample_rest}"
          else
            printf '%s\n' "${staged}"
          fi

        done > sorted_crams.list

        ~{glimpse_phase}  \
            --bam-list sorted_crams.list \
            --reference ~{reference_chunk} \
            --fasta "~{ref_fasta}" \
            --output phase_output.bcf \
            --threads 1

        bcftools index -f phase_output.bcf

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
        File imputed_bcf = "phase_output.bcf"
        File imputed_bcf_index = "phase_output.bcf.csi"
    }

}
