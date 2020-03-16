BIN=${1:-python3}
# $BIN -m spiper get_changed_files pipeline_rnaseq_hisat2_stringtie@file://$PWD TOPLEVEL:main _temp_build/root --plain
set -e
PACKAGE=pipeline_rnaseq_hisat2_stringtie@file://$PWD
$BIN -m spiper get_changed_files $PACKAGE TOPLEVEL:test_job _temp_build/root --plain
$BIN -m spiper run               $PACKAGE TOPLEVEL:test_job _temp_build/root

