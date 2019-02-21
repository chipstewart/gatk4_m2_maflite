task gatk4_m2_maflite_task_1 {
    #Inputs and constants defined here
    String pair_id
    String tumor_id
    String normal_id
    File M2_vcf_file
    String boot_disk_gb = "10"
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions


    command {
        set -euo pipefail
        set -x
        julia --version
        ls -la /opt/src
        ls -la
        ls -la ${M2_vcf_file}
        
        /opt/src/gatk4_m2_maflite.sh ${tumor_id} ${normal_id} ${pair_id} ${M2_vcf_file}

        tar cvfz ${pair_id}.m2_maflite.tar.gz tmp1.tsv ${pair_id}.raw.tsv ${pair_id}.m2.all.maflite.tsv ${pair_id}.m2.pass.maflite.tsv

    }

    output {
        File m2_pass_maflite="${pair_id}.m2.pass.maflite.tsv"
        File m2_all_maflite="${pair_id}.m2.all.maflite.tsv"
        File m2_maflite_tarball="${pair_id}.m2_maflite.tar.gz"
    }

    runtime {
        docker : "chipstewart/gatk4_m2_maflite_task_1:1"
        memory: "${if defined(ram_gb) then ram_gb else '3'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '50'} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }
}

workflow gatk4_m2_maflite {
    call gatk4_m2_maflite_task_1 
}
