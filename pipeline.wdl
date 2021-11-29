version 1.0


task SamtoolsCollate{
    input{
        File bamfile
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
    }
}


task SplitSam{
    input{
        File bamfile
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
    }
}

task ExtractReads{
    input{
        File bamfile
        String sample_id
        Float perc_threshold = 0.85
        Int mapping_quality = 30
        Int telomere_length_threshold = 36
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
        walltime: "4:00"
    }
}

task MergeCsv{
    input{
        Array[File] inputs
        String filename_prefix
    }
    command<<<
        telomere_utils merge_files \
        --inputs ~{sep=" "inputs} --output ~{filename_prefix}.csv.gz
    >>>
    output{
        File out_csv = "~{filename_prefix}.csv.gz"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "4:00"
    }
}


task GetOverlap{
    input{
        File normal_bam
        File normal_csv
        File tumour_csv
        Int binsize = 1000
    }
    command<<<
        telomere_utils get_overlap \
        --normal_bam ~{normal_bam} --normal_data ~{normal_csv} \
        --tumour_data ~{tumour_csv} --output overlapping.csv.gz \
        --bin_counts bin_counts.csv.gz --binsize ~{binsize}
    >>>
    output{
        File outfile = "overlapping.csv.gz"
        File bin_counts = "bin_counts.csv.gz"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "4:00"
    }
}




workflow TelomereWorkflow{
    input{
        File normal_bam
        File tumour_bam
        String normal_sample_id
        String tumour_sample_id
        Float perc_threshold = 0.85
        Int mapping_quality = 30
        Int telomere_length_threshold = 36
        Int binsize = 1000
    }

    call SamtoolsCollate as normal_collate{
        input:
            bamfile = normal_bam
    }
    call SamtoolsCollate as tumour_collate{
        input:
            bamfile = tumour_bam
    }
    call SplitSam as normal_split{
        input:
        bamfile = normal_collate.collated_bam
    }
    call SplitSam as tumour_split{
        input:
        bamfile = tumour_collate.collated_bam
    }

    scatter (bamfile in  normal_split.bamfiles){
        call ExtractReads as extract_normal{
            input:
                bamfile = bamfile,
                sample_id = normal_sample_id,
                perc_threshold=perc_threshold,
                mapping_quality=mapping_quality,
                telomere_length_threshold=telomere_length_threshold
        }
    }

    scatter (bamfile in  tumour_split.bamfiles){
        call ExtractReads as extract_tumour{
            input:
                bamfile = bamfile,
                sample_id = tumour_sample_id,
                perc_threshold=perc_threshold,
                mapping_quality=mapping_quality,
                telomere_length_threshold=telomere_length_threshold
        }
    }

    call MergeCsv as merge_normal{
        input:
            inputs = extract_normal.outcsv,
            filename_prefix = "normal_telomeres"
    }
    call MergeCsv as merge_tumour{
        input:
            inputs = extract_tumour.outcsv,
            filename_prefix = "tumour_telomeress"
    }

    call GetOverlap as overlap{
        input:
            normal_bam = normal_bam,
            normal_csv = merge_normal.out_csv,
            tumour_csv = merge_tumour.out_csv,
            binsize=binsize
    }

    output{
        File tumour_data = merge_tumour.out_csv
        File normal_data = merge_normal.out_csv
        File overlap_data = overlap.outfile
        File bin_counts = overlap.bin_counts
    }

}
