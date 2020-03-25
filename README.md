[![Build Status](https://travis-ci.com/shouldsee/pipeline-rnaseq-hisat2-stringtie.svg?branch=master)](https://travis-ci.com/shouldsee/pipeline-rnaseq-hisat2-stringtie)

# pipeline-rnaseq-hisat2-stringtie

```
Comment: This is the default pair-end RNA-Seq pipeline with hisat2 during Feng's time at SLCU.
Author: Feng Geng
Last Modified: 2020/03/16
version: 0.0.6
Example Usage: 

	BIN=python3
	PACKAGE=pipeline_rnaseq_hisat2_stringtie@https://github.com/shouldsee/pipeline-rnaseq-hisat2-stringtie/tarball/master

	## install dependencies
	$BIN -m pip install spiper@https://github.com/shouldsee/spiper/tarball/master --user --upgrade
	curl -sL https://raw.githubusercontent.com/shouldsee/spiper/master/scripts/install_singular.sh | bash -s $HOME/.local

	## test example runs
	$BIN -m spiper get_changed_files $PACKAGE TOPLEVEL:test_job _temp_build/root --plain
	$BIN -m spiper run               $PACKAGE TOPLEVEL:test_job _temp_build/root


	time $BIN -m spiper run $PACKAGE TOPLEVEL:workflow \
		$HOME/run_temp/0306.201RS6/root \
		$HOME/_hisat2_cache/Bd21_v3-1 \
		$HOME/ref/Bd21_v3-1/assembly/Bdistachyon_314_v3.0.fa \
		$HOME/ref/Bd21_v3-1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf \
		$HOME/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/R1_001.fastq.gz \
		$HOME/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/R2_001.fastq.gz \
		8

	ls $HOME/run_temp/0306.201RS6/ -lhtr

Example output:

	total 15G
	drwxrwxr-x  4 feng feng 4.0K Mar 16 15:51 root.job_trimmomatic.singularity_temp
	-rw-rw-r--  1 feng feng 4.0G Mar 16 15:56 root.job_trimmomatic.fastq1
	-rw-rw-r--  1 feng feng 3.9G Mar 16 15:56 root.job_trimmomatic.fastq2
	-rw-rw-r--  1 feng feng 5.7M Mar 16 15:56 root.job_trimmomatic.fastq1.fail
	-rw-rw-r--  1 feng feng 1.1M Mar 16 15:56 root.job_trimmomatic.fastq2.fail
	-rw-rw-r--  1 feng feng 1.7K Mar 16 15:56 root.job_trimmomatic.log
	-rw-rw-r--  1 feng feng 5.6K Mar 16 15:56 root.job_trimmomatic.cmd
	drwxrwxr-x  4 feng feng 4.0K Mar 16 15:56 root.job_hisat2_align.singularity_temp
	-rw-rw-r--  1 feng feng 4.4M Mar 16 16:06 root.job_hisat2_align.log
	-rw-rw-r--  1 feng feng 2.8G Mar 16 16:06 root.job_hisat2_align.bam.unsorted
	-rw-rw-r--  1 feng feng 376K Mar 16 16:08 root.job_hisat2_align.log.2
	-rw-rw-r--  1 feng feng 1.9G Mar 16 16:09 root.job_hisat2_align.bam
	drwxrwxr-x  2 feng feng 4.0K Mar 16 16:09 root.job_hisat2_align.bam.sort_temp
	-rw-rw-r--  1 feng feng  11K Mar 16 16:09 root.job_hisat2_align.cmd
	drwxrwxr-x  4 feng feng 4.0K Mar 16 16:09 root.job_picard_dedup.singularity_temp
	-rw-rw-r--  1 feng feng 1.6G Mar 16 16:26 root.job_picard_dedup.bam
	-rw-rw-r--  1 feng feng 6.5K Mar 16 16:26 root.job_picard_dedup.log
	-rw-rw-r--  1 feng feng 383K Mar 16 16:26 root.job_picard_dedup.bam.bai
	-rw-rw-r--  1 feng feng  25K Mar 16 16:26 root.job_picard_dedup.cmd_log
	drwxrwxr-x  4 feng feng 4.0K Mar 16 16:26 root.job_stringtie_count.singularity_temp
	-rw-rw-r--  1 feng feng 1.8M Mar 16 16:33 root.job_stringtie_count.count
	-rw-rw-r--  1 feng feng  17M Mar 16 16:33 root.job_stringtie_count.cmd
	drwxrwxr-x 37 feng feng 4.0K Mar 16 16:33 _spiper
```

