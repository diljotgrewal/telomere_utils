version 1.0


task SamtoolsCollate{
    input{
        File bamfile
        File baifile
        String? singularity_image
    }
    command<<<
        mkdir tempdir
        samtools collate -u ~{bamfile} tempdir/output
        samtools view -h tempdir/output.bam -o tempdir/output.sam
    >>>
    output{
        File collated_bam = "tempdir/output.sam"
    }
    runtime{
        memory: "25 GB"
        cpu: 1
        walltime: "24:00"
        singularity: '~{singularity_image}'
        docker: 'quay.io/wgspipelinetest/telomere_utils:v0.0.1'
    }
}


task SplitSam{
    input{
        File bamfile
        String? singularity_image
    }
    command<<<
        telomere_utils split_bam --infile ~{bamfile} --outdir tempdir/
    >>>
    output{
        Array[File] bamfiles = glob("tempdir/*.sam")
    }
    runtime{
        memory: "6 GB"
        cpu: 1
        walltime: "10:00"
        singularity: '~{singularity_image}'
        docker: 'quay.io/wgspipelinetest/telomere_utils:v0.0.1'
    }
}

task ExtractReads{
    input{
        File bamfile
        String sample_id
        Float perc_threshold = 0.85
        Int mapping_quality = 30
        Int telomere_length_threshold = 36
        String? singularity_image
    }
    command<<<
        telomere_utils extract_reads --input ~{bamfile} \
        --outbam output.bam --outcsv output.csv \
        --sample_id ~{sample_id} --perc_threshold ~{perc_threshold} \
        --mapping_quality ~{mapping_quality} \
        --telomere_length_threshold ~{telomere_length_threshold}
    >>>
    output{
        File outbam = "output.bam"
        File outcsv = "output.csv"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "24:00"
        singularity: '~{singularity_image}'
        docker: 'quay.io/wgspipelinetest/telomere_utils:v0.0.1'
    }
}

task MergeCsv{
    input{
        Array[File] inputs
        String filename_prefix
        String? singularity_image
    }
    command<<<
        telomere_utils merge_files \
        --inputs ~{sep=" "inputs} --output ~{filename_prefix}.csv.gz
    >>>
    output{
        File out_csv = "~{filename_prefix}.csv.gz"
    }
    runtime{
        memory: "60 GB"
        cpu: 1
        walltime: "4:00"
        singularity: '~{singularity_image}'
        docker: 'quay.io/wgspipelinetest/telomere_utils:v0.0.1'
    }
}


task GetOverlap{
    input{
        File normal_bam
        File normal_bai
        File normal_csv
        File tumour_csv
        String chromosome
        Int binsize = 1000
        String? singularity_image
    }
    command<<<
        telomere_utils get_overlap \
        --normal_bam ~{normal_bam} --normal_data ~{normal_csv} \
        --tumour_data ~{tumour_csv} --output overlapping.csv.gz \
        --bin_counts bin_counts.csv.gz --binsize ~{binsize} \
        --chromosome ~{chromosome}
    >>>
    output{
        File outfile = "overlapping.csv.gz"
        File bin_counts = "bin_counts.csv.gz"
    }
    runtime{
        memory: "60 GB"
        cpu: 1
        walltime: "24:00"
        singularity: '~{singularity_image}'
        docker: 'quay.io/wgspipelinetest/telomere_utils:v0.0.1'
    }
}




workflow TelomereWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        String normal_sample_id
        String tumour_sample_id
        Array[String] chromosomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
        Float tumour_perc_threshold = 0.85
        Int tumour_mapping_quality = 30
        Int tumour_telomere_length_threshold = 36
        Float normal_perc_threshold = 0.85
        Int normal_mapping_quality = 30
        Int normal_telomere_length_threshold = 36
        Int binsize = 1000
        String? singularity_image
    }

    call SamtoolsCollate as normal_collate{
        input:
            bamfile = normal_bam,
            baifile = normal_bai,
            singularity_image = singularity_image
    }
    call SamtoolsCollate as tumour_collate{
        input:
            bamfile = tumour_bam,
            baifile = tumour_bai,
            singularity_image = singularity_image
    }
    call SplitSam as normal_split{
        input:
            bamfile = normal_collate.collated_bam,
            singularity_image = singularity_image
    }
    call SplitSam as tumour_split{
        input:
            bamfile = tumour_collate.collated_bam,
            singularity_image = singularity_image
    }

    scatter (bamfile in  normal_split.bamfiles){
        call ExtractReads as extract_normal{
            input:
                bamfile = bamfile,
                sample_id = normal_sample_id,
                perc_threshold=normal_perc_threshold,
                mapping_quality=normal_mapping_quality,
                telomere_length_threshold=normal_telomere_length_threshold,
                singularity_image = singularity_image
        }
    }

    scatter (bamfile in  tumour_split.bamfiles){
        call ExtractReads as extract_tumour{
            input:
                bamfile = bamfile,
                sample_id = tumour_sample_id,
                perc_threshold=tumour_perc_threshold,
                mapping_quality=tumour_mapping_quality,
                telomere_length_threshold=tumour_telomere_length_threshold,
                singularity_image = singularity_image
        }
    }

    call MergeCsv as merge_normal{
        input:
            inputs = extract_normal.outcsv,
            filename_prefix = "normal_telomeres",
            singularity_image = singularity_image
    }
    call MergeCsv as merge_tumour{
        input:
            inputs = extract_tumour.outcsv,
            filename_prefix = "tumour_telomeres",
            singularity_image = singularity_image
    }

    scatter(chrom in chromosomes){
        call GetOverlap as overlap{
            input:
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                normal_csv = merge_normal.out_csv,
                tumour_csv = merge_tumour.out_csv,
                binsize=binsize,
                chromosome = chrom,
                singularity_image = singularity_image
        }
    }

    call MergeCsv as merge_overlaps{
        input:
            inputs = overlap.outfile,
            filename_prefix = "telomeres_overlap",
            singularity_image = singularity_image
    }

    call MergeCsv as merge_bin_counts{
        input:
            inputs = overlap.bin_counts,
            filename_prefix = "tumour_telomeres_bin_counts",
            singularity_image = singularity_image
    }

    output{
        File tumour_data = merge_tumour.out_csv
        File normal_data = merge_normal.out_csv
        File overlap_data = merge_overlaps.out_csv
        File bin_counts = merge_bin_counts.out_csv
    }

}