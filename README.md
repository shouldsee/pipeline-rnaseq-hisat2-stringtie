[![Build Status](https://travis-ci.com/shouldsee/pipeline-rnaseq-hisat2-stringtie.svg?branch=master)](https://travis-ci.com/shouldsee/pipeline-rnaseq-hisat2-stringtie)

# pipeline-rnaseq-hisat2-stringtie

```
Comment: This is the default pair-end RNA-Seq pipeline with hisat2 during Feng's time at SLCU.
Author: Feng Geng
Last Modified: 2020/03/16
Example Usage: 

	BIN=python3
	PACKAGE=pipeline_rnaseq_hisat2_stringtie@https://github.com/shouldsee/pipeline-rnaseq-hisat2-stringtie/tarball/master

	## install dependencies
	$BIN -m pip install spiper@https://github.com/shouldsee/spiper/tarball/master --user --upgrade
	curl -sL https://raw.githubusercontent.com/shouldsee/spiper/master/scripts/install_singular.sh | bash -s $HOME/.local

	## test example runs
	$BIN -m spiper get_changed_files $PACKAGE TOPLEVEL:test_job _temp_build/root --plain
	$BIN -m spiper run               $PACKAGE TOPLEVEL:test_job _temp_build/root

	$BIN -m spiper run $PACKAGE TOPLEVEL:workflow \
		$HOME/run_temp/0306.201RS6/root \
		$HOME/_hisat2_cache/Bd21_v3-1 \
		$HOME/ref/Bd21_v3-1/assembly/Bdistachyon_314_v3.0.fa \
		$HOME/ref/Bd21_v3-1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf \
		$HOME/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/R1_001.fastq.gz \
		$HOME/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/R2_001.fastq.gz \
		8
```

