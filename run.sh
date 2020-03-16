BIN=${1:-python3}
# $BIN -m spiper get_changed_files pipeline_rnaseq_hisat2_stringtie@file://$PWD TOPLEVEL:main _temp_build/root --plain
set -e
$BIN -m spiper get_changed_files $PWD pipeline_rnaseq_hisat2_stringtie:main _temp_build/root --plain
$BIN -m spiper run $PWD pipeline_rnaseq_hisat2_stringtie:main _temp_build/root

