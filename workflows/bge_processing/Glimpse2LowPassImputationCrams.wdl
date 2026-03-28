version 1.0

workflow Glimpse2LowPassImputationCrams {

    input {

        String pipeline_version = "0.0.1"

        Array[File] crams
        Array[File] crais
        Array[String] sample_ids
        File reference_chunks
        String output_basename

        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict

        String docker

    }

    scatter (reference_chunk in read_lines(reference_chunks)) {

        call GlimpsePhase {
            input:
                crams = crams,
                crais = crais,
                sample_ids = sample_ids,
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
             output_basename = output_basename,
             docker = docker
    }

    output {
        File imputed_vcf = GlimpseLigate.imputed_vcf
        File imputed_vcf_index = GlimpseLigate.imputed_vcf_index
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

        String glimpse_ligate= "/usr/local/bin/GLIMPSE2_ligate"
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
 
task GlimpsePhase{

    input {

        File reference_chunk

        Array[File] crams
        Array[File] crais
        Array[String] sample_ids

        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict

        Int mem_gb = 16 
        Int disk_size_gb = 100
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 2

        String docker
        String glimpse_phase= "/usr/local/bin/GLIMPSE2_phase"

    }

    command <<<

        ~{glimpse_phase}  \
            --bam-list ~{write_lines(crams)} \
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

