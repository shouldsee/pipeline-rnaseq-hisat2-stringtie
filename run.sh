BIN=${1:-python3}
# $BIN -m spiper get_changed_files pipeline_rnaseq_hisat2_stringtie@file://$PWD TOPLEVEL:main _temp_build/root --plain
set -e
PACKAGE=pipeline_rnaseq_hisat2_stringtie@file://$PWD
$BIN -m spiper get_changed_files $PACKAGE TOPLEVEL:test_job _temp_build/root --plain
$BIN -m spiper run               $PACKAGE TOPLEVEL:test_job _temp_build/root

exit 0

PACKAGE=pipeline_rnaseq_hisat2_stringtie@https://github.com/shouldsee/pipeline-rnaseq-hisat2-stringtie/tarball/master
time python3 -m spiper run $PACKAGE TOPLEVEL:workflow \
	./201RS6.0306/root \
	~/_hisat2_cache/Bd21_v3-1 \
	/home/feng/ref/Brachypodium_Bd21_v3-1/genome.fa \
	/home/feng/ref/Brachypodium_Bd21_v3-1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf \
	~/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/*R1*.fastq.gz \
	~/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/*R2*.fastq.gz \
	4

python3 -m spiper get_changed_files $PACKAGE TOPLEVEL:workflow \
	./201RS6.0306/root \
	~/_hisat2_cache/Bd21_v3-1 \
	/home/feng/ref/Bd21_v3-1/assembly/Bdistachyon_314_v3.0.fa \
	/home/feng/ref/Bd21_v3-1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf \
	~/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/*R1*.fastq.gz \
	~/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/*R2*.fastq.gz \
	4

python3 -m spiper run $PACKAGE TOPLEVEL:workflow \
	./201RS6.0306/root \
	~/_hisat2_cache/Bd21_v3-1 \
	/home/feng/ref/Bd21_v3-1/assembly/Bdistachyon_314_v3.0.fa \
	/home/feng/ref/Bd21_v3-1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf \
	~/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/*R1*.fastq.gz \
	~/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/*R2*.fastq.gz \
	4



python3 -m spiper get_changed_files $PACKAGE TOPLEVEL:workflow ./201RS6.0306/root ~/_hisat2_cache/Bd21_v3-1 /home/feng/ref/Brachypodium_Bd21_v3-1/genome.fa /home/feng/ref/Brachypodium_Bd21_v3-1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf ~/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/*R1*.fastq.gz ~/envs/upGeo/brachy-meta/rawfile.workdir/201RS6/combined_valid_fastq/*R2*.fastq.gz 4