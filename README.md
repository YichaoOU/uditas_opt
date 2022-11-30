# uditas_opt
Optimized version of the original UDITAS program: https://github.com/editasmedicine/uditas

# CHANGES

1. HBG promoters are almost identical, reads can be mapped to both `wt` or other amplicons like `large_deletion`. The change was made to prioritize wt reads by post-processing bowtie2 bam file; any reads mapped to wt reads will not be mapped to other amplicons. see function `fix_multi_mapped_wt_reads_bam`

2. New function `create_amplicon_user` return `amplicon_list`, derived from `create_amplicon`, so that it accepts a list of break sites to quantify translocation. Modified `create_amplicon` to take the list `amplicon_list` and create amplicon sequences.

3. Modified `align_genome_local` to run SV prediction using `Socrates`

4. Modified `cutadapt` parameters `-m 20`.

5. `sample_info.csv`, each sample needs to have 2 gRNA locations.


# Note to usage

1. Input fastq file names are fixed

```
	r1_fastq = os.path.join(dir_sample, 'Undetermined_S0_R1_001.fastq.gz')
	r2_fastq = os.path.join(dir_sample, 'Undetermined_S0_R2_001.fastq.gz')
	i1_fastq = os.path.join(dir_sample, 'Undetermined_S0_I1_001.fastq.gz')
	i2_fastq = os.path.join(dir_sample, 'Undetermined_S0_I2_001.fastq.gz')
```

2. Input sample info file is fixed

```

sample_info_filename = os.path.join(dir_sample, 'sample_info.csv')
```

# How this pipeline work:

1. Demultiplexing, not parallelized, slow, mismatch=1. Use `-skip_demultiplexing 1` to skip

2. For each sample in `sample_info.csv`, do the following

3. trim fastq. Use `-skip_trimming 1` to skip

4. Perform read alignment to the genome, bowtie2 local, then identify SVs using Socrates (`align_genome_local`, if `socrates_path` is not None, SV will be performed). Next, any SVs overlapping with the given gRNA location (`parse_socrates_output`) will be used to create amplicon sequence (`create_amplicon_user`). 

5. Create amplicon sequences to quantify editing outcomes.  `create_amplicon` -> `align_amplicon` -> `analyze_alignments` & `analyze_alignments_all_amplicons`, the latter is just used to get a total number of mapped reads.

6. For unmapped reads,  `extract_unmapped_reads_amplicons` -> `analyze_alignments_genome_global` to identify mis-priming. Not sure why do they use global mapping.

7. Summary of everything is in this file `results_summary_pivot.xlsx`
