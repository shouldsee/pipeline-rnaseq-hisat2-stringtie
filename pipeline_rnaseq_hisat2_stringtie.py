__doc__ = '''
Comment: This is the default pair-end RNA-Seq pipeline with hisat2 during Feng's time at SLCU.
Author: Feng Geng
Last Modified: 2020/03/18
Version: 0.0.5
Example Usage: 

	BIN=python3
	PACKAGE=pipeline_rnaseq_hisat2_stringtie@https://github.com/shouldsee/pipeline-rnaseq-hisat2-stringtie/tarball/master

	## install dependencies
	$BIN -m pip install spiper@https://github.com/shouldsee/spiper/tarball/master --user --upgrade
	curl -sL https://raw.githubusercontent.com/shouldsee/spiper/master/scripts/install_singular.sh | bash -s $HOME/.local

	## test example runs
	## (test_run)
	$BIN -m spiper get_all_files $PACKAGE TOPLEVEL:test_job _temp_build/root
	$BIN -m spiper get_changed_files $PACKAGE TOPLEVEL:test_job _temp_build/root --plain
	$BIN -m spiper get_all_deps      $PACKAGE TOPLEVEL:test_job _temp_build/root
	$BIN -m spiper run               $PACKAGE TOPLEVEL:test_job _temp_build/root

	## (production_run)
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
	ls: (production_run)
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

	get_all_files: (test_run)
		[File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.get_fasta.fasta'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.get_fasta.cmd'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.get_genepred.genepred'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.get_genepred.gtf'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.get_genepred.cmd'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.hisat2/wuhan-ncov19.job_hisat2_index.log'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.hisat2/wuhan-ncov19.job_hisat2_index.cmd'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_trimmomatic.fastq1'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_trimmomatic.fastq2'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_trimmomatic.log'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_trimmomatic.cmd'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_hisat2_align.bam'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_hisat2_align.log'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_hisat2_align.cmd'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_picard_dedup.bam'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_picard_dedup.log'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_picard_dedup.cmd_log'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_bam2bw-picard_dedup_bam.bw'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_bam2bw-picard_dedup_bam.cmd'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_stringtie_count.count'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/_temp_build/root.sample1.job_stringtie_count.cmd')
		 ]

	get_all_deps:  (test_run) (spiper>=0.1.1)
		[Depend('bin://curl'),
		 Depend('bin://gzip'),
		 Depend('docker://quay.io/biocontainers/hisat2:2.1.0--py36hc9558a2_4'),
		 Depend('docker://quay.io/biocontainers/picard:2.21.9--0'),
		 Depend('docker://quay.io/biocontainers/samtools:1.10--h9402c20_2'),
		 Depend('docker://quay.io/biocontainers/stringtie:2.1.1--hc900ff6_0'),
		 Depend('docker://quay.io/biocontainers/trimmomatic:0.35--6'),
		 Depend('docker://quay.io/biocontainers/ucsc-genepredtogtf:377--h35c10e6_2'),
		 Depend('docker://quay.io/wtsicgp/cgpbigwig:1.1.0'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/test_data/test_R1_.fastq'),
		 File('/home/user/repos/pipeline-rnaseq-hisat2-stringtie/test_data/test_R2_.fastq'),
		 HttpResponseContentHeader({"method"="head","url"="https://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Wuhan_seafood_market_pneumonia_virus/bigZips/chromFa.tar.gz"}),
		 HttpResponseContentHeader({"method"="head","url"="https://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Wuhan_seafood_market_pneumonia_virus/database/ncbiGene.txt.gz"})]    

'''

from spiper.types import File,TempFile, Prefix, Default
from spiper.types import job_result
from spiper.types import Depend
from spiper.runner import list_flatten_strict, job_from_func

from spiper.shell import LoggedSingularityCommandList, LoggedShellCommand, LoggedSingularityCommand
from path import Path
import spiper
assert spiper.VERSION >= '0.1.1',spiper.VERSION

from spiper.types import Concat
from spiper.types import Flow
from spiper.types import resolve_spiper
from spiper.types import LoggedShellCommand
from spiper.types import CopyFile
from spiper.types import RPO


def job_trimmomatic(
	self, prefix,
	FASTQ_FILE_1 = File, 
	FASTQ_FILE_2 = File, 
	THREADS_ = int,
	_IMAGE = Depend('docker://quay.io/biocontainers/trimmomatic:0.35--6'),
	_output = [
		File('fastq1'),
		File('fastq2'),
		File('log'),
		File('cmd'),
		],
	):    
		_ = '''
		trimmomatic PE -threads 4 -phred33 
		/home/feng/temp/187R/187R-S1-2018_06_27_14:02:08/809_S1_R1_raw.fastq 
		/home/feng/temp
	/187R/187R-S1-2018_06_27_14:02:08/809_S1_R2_raw.fastq 
	809_S1_R1_raw_pass.fastq 
	809_S1_R1_raw_fail.fastq 
	809_S1_R2_raw_pass.fastq 
	809_S1_R2_raw_fail.fastq 
	ILLUMINACLIP:/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa
	:6:30:10 LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15
		'''
		# _out = get_output_files(self, prefix, _output)

		CMD = [
		'trimmomatic',
		'PE','-threads', str(THREADS_//2), 
		'-phred33',
		File( FASTQ_FILE_1 ),
		File( FASTQ_FILE_2 ),
		File( self.output['fastq1'] ),
		File( self.output['fastq1'] + '.fail'),
		File( self.output['fastq2'] ),
		File( self.output['fastq2'] + '.fail'),
		'ILLUMINACLIP:'
		'/usr/local/share/trimmomatic-0.35-6/adapters/TruSeq3-PE-2.fa'
		':6:30:10',
		'LEADING:3',
		'TRAILING:3',
		'MINLEN:36',
		'SLIDINGWINDOW:4:15',
		'&>', 
		File( self.output['log'])
		]
		res = LoggedSingularityCommand(self.prefix_named, CMD, _IMAGE, self.output['cmd'])
		return self
		# return job_result( None, CMD, self.output)

def job_hisat2_index(
	self,prefix, 
	FASTA_FILE = File,
	THREADS_  = int,
	_IMAGE    = Depend("docker://quay.io/biocontainers/hisat2:2.1.0--py36hc9558a2_4"),
	_output   = [
		Prefix('index_prefix'), 
		File('log'),
		File('cmd'),
	],
	):

	CMD = [
	'hisat2-build',
	'-p',str(THREADS_),
	 File(  FASTA_FILE),
	 Prefix(self.output.index_prefix),
	 '&>', 
	 File(  self.output.log),
	 ]
	res = LoggedSingularityCommand(self.prefix_named, CMD, _IMAGE, self.output.cmd)
	return self


def job_hisat2_align(
	self,prefix,
	INDEX_PREFIX = Prefix,
	FASTQ_FILE_1 = File,
	FASTQ_FILE_2 = File,
	hisat2_args = list,
	THREADS_ = int,
	_IMAGE   = Depend("docker://quay.io/biocontainers/hisat2:2.1.0--py36hc9558a2_4"),
	_IMAGE_SAMTOOLS = Depend("docker://quay.io/biocontainers/samtools:1.10--h9402c20_2"),
	_output = [
		File('bam'),
		File('log'),
		File('cmd'),
	]
	):
	# _out = get_output_files(self,prefix,_output)
	results = []
	cmd1 = CMD = [
	 'hisat2',
	 # hisat2_args,
	 '-x', Prefix(INDEX_PREFIX),
	 '-1', File( FASTQ_FILE_1),
	 '-2', File( FASTQ_FILE_2),
	 # '-U', File( FASTQ_FILE_1),
	 # ['-2',File( FASTQ_FILE_2) ] if FASTQ_FILE_2 else [],
	 '-S', '/dev/stdout',
	 '--threads', str( max(1, THREADS_-1) ),
	 hisat2_args or [
	 '--no-mixed',
	 '--rna-strandness','RF',
	 '--dta',
	 '--fr'],
	 '2>', File( self.output.log),
	]
	'''
	singularity --verbose --debug exec docker://python:2.7.17-alpine python -V
	singularity shell docker://python:2.7.17-alpine python -V
	'''
	# res = LoggedSingularityCommand(CMD, _IMAGE, self.output.cmd)

	# results.append(job_result( None, CMD, self.output))
	# _ = '''
	# samtools view /home/feng/temp/187R/187R-S1-2018_06_27_14:02:08/809_S1.sam -b --threads 4 -o 809_S1.bam
	# '''
	cmd2 = CMD = [    
		'samtools','view','-bS','/dev/stdin',
		'--threads',str( 1 ),
		'-o', 
		( self.output.bam+'.unsorted'),
	]
	# res = LoggedSingularityCommand(CMD, _IMAGE_SAMTOOLS, self.output.cmd)

	cmd3 = CMD = [
		'samtools','sort', ( self.output.bam + '.unsorted'),
		'--threads', str( THREADS_ ),
		'-o', ( self.output.bam),
		'-T', File(self.output.bam+'.sort_temp/').makedirs_p().check_writable(),
	]

	CMD = [
		# 'PIPE=$(mktemp -u);mkfifo $PIPE;exec 3<>$PIPE ;rm $PIPE;',
		LoggedSingularityCommandList(self.prefix_named, cmd1, _IMAGE,),'|',
		LoggedSingularityCommandList(self.prefix_named, cmd2, _IMAGE_SAMTOOLS),'&&',
		LoggedSingularityCommandList(self.prefix_named, cmd3, _IMAGE_SAMTOOLS),
		 # extra_files = [File(self.output.bam.dirname())]),
		# LoggedSingularityCommandList(cmd3, _IMAGE_SAMTOOLS, extra_files = [File(self.output.bam.dirname())]),
		# LoggedSingularityCommandList([cmd3,'&&','df',File(self.output.bam.dirname())], _IMAGE_SAMTOOLS, 
		#     extra_files = [File(self.output.bam.dirname())]),
	]
	res = LoggedShellCommand(CMD, self.output.cmd)
	# (self.output.bam+'.sam').unlink_p()
	# (self.output.bam+'.unsorted').unlink_p()

	# res = LoggedSingularityCommand(CMD, _IMAGE_SAMTOOLS, self.output.cmd)
	return self


def job_stringtie_count(self, prefix,
	BAM_FILE = File,
	GTF_FILE = File,
	THREADS_ = int,    
	_IMAGE = Depend('docker://quay.io/biocontainers/stringtie:2.1.1--hc900ff6_0'),
	_output = ['count','cmd']
	):
	_= '''
	Example run:
		stringtie 
		-p 4 
		--rf 809_S1.bam 
		-G /home/feng/ref/Arabidopsis_thaliana_TAIR10/annotation/genes.gtf 
		-o 809_S1.stringtie.gtf 
		-A 809_S1.stringtie.count &> 809_S1.stringtie.log
	'''
	CMD = [
	'stringtie',
	'-p', str(THREADS_), File(BAM_FILE),
	'--rf',
	'-G', File(GTF_FILE),
	'-A', File(self.output.count),
	]
	res = LoggedSingularityCommand(self.prefix_named, CMD, _IMAGE, self.output.cmd)

def job_picard_dedup(
	self,prefix,
	bam_file = File,
	THREADS_ = int,
	_IMAGE = Depend('docker://quay.io/biocontainers/picard:2.21.9--0'),
	_IMAGE_SAMTOOLS = Depend("docker://quay.io/biocontainers/samtools:1.10--h9402c20_2"),
	_output= ['bam','log','cmd_log'],):
		CMD = [
		'picard',
		'MarkDuplicates',
		Concat('I=',File(bam_file)),
		Concat('O=',File(self.output.bam)),
		Concat('M=',File(self.output.log)),
		# Concat('TMP_DIR=',File(self.output.bam+'.picard_temp').makedirs_p().check_writable()),
		'REMOVE_DUPLICATES=true',
		]
		res = LoggedSingularityCommand(self.prefix_named, CMD, _IMAGE, self.output.cmd_log,)
		res = LoggedSingularityCommand(
			self.prefix_named,
			# prefix,
			['samtools','index',self.output.bam],
			_IMAGE_SAMTOOLS, 
			self.output.cmd_log,mode='a',
			extra_files = [self.output.bam+'.bai'])

import json,math
def job_bam2bw_cpm(self,prefix, 
	bam_file= File,
	bam_qc_file = File,
	THREADS_ = int,
	_image = Depend('docker://quay.io/shouldsee/cgpbigwig:b024993'),
	# _image = Depend('docker://quay.io/wtsicgp/cgpbigwig:1.1.0'),
	_output=['cpm_bw','cmd'],
	):
	'''
	#### set scale_log10==0. to disable rescaling
	'''
	assert (bam_file+'.bai').isfile()

	scale_log10 = math.log10(1.E6 / max(1,
			json.loads(open(bam_qc_file,'r').read())['counts.uniq_mapped.sum']
		))
	CMD = ['bam2bw','-S',str(scale_log10),'-i',bam_file, '-o', self.output.cpm_bw]
	LoggedSingularityCommand(self.prefix_named, CMD, _image, self.output.cmd,extra_files = [bam_file+'.bai'])




# DATA_DICT['counts']['FWD'] = int(shellexec('samtools view -c -F 0x10 -F 0x100 -F 0x4 ALIGNMENT-sorted.bam').strip())
# DATA_DICT['counts']['REV'] = int(shellexec('samtools view -c -f 0x10 -F 0x100 -F 0x4 ALIGNMENT-sorted.bam').strip())
# DATA_DICT['counts']['SUM'] = DATA_DICT['counts']['FWD'] + DATA_DICT['counts']['REV']
# DATA_DICT['texts'] = pyext._DICT_CLASS()
# DATA_DICT['texts']['flagstat'] = shellexec('samtools flagstat ALIGNMENT-sorted.bam')
# DATA_DICT['texts']['quickcheck'] = shellexec('samtools quickcheck ALIGNMENT-sorted.bam 2>&1')

'''
# example flagstat
36143496 + 0 in total (QC-passed reads + QC-failed reads)
2174156 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
30426164 + 0 mapped (84.18% : N/A)
33969340 + 0 paired in sequencing
16984670 + 0 read1
16984670 + 0 read2
27663462 + 0 properly paired (81.44% : N/A)
28252008 + 0 with itself and mate mapped ( -F0x100 -F 0x4 )
0 + 0 singletons (0.00% : N/A)
437508 + 0 with mate mapped to a different chr
437508 + 0 with mate mapped to a different chr (mapQ>=5)
'''
import collections,json
def job_bam_qc(self,prefix,
	bam_file = File,
	THREADS_ = int,
	_image = Depend("docker://quay.io/biocontainers/samtools:1.10--h9402c20_2"),
	_output = ['cmd','data_json']):
	DATA_DICT = collections.OrderedDict()
	# DATA_DICT['counts'] = collections.OrderedDict()
	cmd_runned , stdout = (LoggedSingularityCommand(self.prefix_named, [
		'bash -euc "{ ',
		'samtools view -c -f 0x4 ',bam_file,';',  ## UNMAPPED
		'samtools view -c -F0x10 -F0x100 -F0x4 ',bam_file,';', ### FWD_UNIQ_MAPPED
		'samtools view -c -f0x10 -F0x100 -F0x4 ',bam_file,';', ### REV_UNIQ_MAPPED
		'}"',
		], _image, self.output.cmd, extra_files = [bam_file+'.bai']))
	sp = stdout.splitlines()
	assert len(sp)==3;
	DATA_DICT['version'] = '0.0.1'
	DATA_DICT['counts.unmapped'] = int(sp[0])
	DATA_DICT['counts.uniq_mapped.fwd'] = int(sp[1])
	DATA_DICT['counts.uniq_mapped.rev'] = int(sp[2])
	DATA_DICT['counts.uniq_mapped.sum'] = int(sp[1]) + int(sp[2])
	DATA_DICT['filename'] = str(bam_file)
	with open(self.output.data_json,'w') as f: json.dump(DATA_DICT,f,indent=2)


@Flow
def workflow(self, prefix, 

	hisat2_cache_prefix = str,
	genome_fasta = File, 
	genome_gtf_file = File,

	fastq1 = File,
	fastq2 = File,

	THREADS_ = int,
	_output=[]
	):
	# print
	# assert 0,repr((hisat2_cache_prefix))
	# self.data = {}
	# self.data['index'] = 
	curr = self.runner(
		job_hisat2_index, 
		hisat2_cache_prefix,
		genome_fasta,
		THREADS_,
		)
	# self.data['trimmed'] = 
	curr = self.runner(
		job_trimmomatic,
		prefix,
		fastq1,
		fastq2,
		THREADS_,
		)
	curr = self.runner(
		job_hisat2_align,
		prefix,
		self.subflow['job_hisat2_index'].output.index_prefix,
		self.subflow['job_trimmomatic'].output.fastq1,
		self.subflow['job_trimmomatic'].output.fastq2,
		[],
		THREADS_,
		)
	self.runner(
		job_picard_dedup,
		prefix,
		self.subflow['job_hisat2_align'].output.bam,
		THREADS_,
		)

	self.config_runner(tag='picard_dedup_bam')(
		job_bam_qc,
		prefix,
		self.subflow['job_picard_dedup'].output.bam,
		THREADS_,
		)
	self.config_runner(tag='picard_dedup_bam')(
		job_bam2bw_cpm,
		prefix,
		self.subflow['job_picard_dedup'].output.bam,
		self.subflow['job_bam_qc-picard_dedup_bam'].output.data_json,
		1,
		)

	self.runner(
		job_stringtie_count,
		prefix,
		self.subflow['job_picard_dedup'].output.bam,
		genome_gtf_file,
		THREADS_,
		)

	assert File(__file__).isfile(),'Cannot find source file using __file__:%r'%__file__
	self.runner(
		CopyFile, 
		self.prefix_named+'.source.py',
		__file__)

	return self


# from spiper.types import Caller, DirtyKey, rgetattr
# import shutil
# def LinkEvent(self, 
#     prefix, input=File, 
#     _single_file = 1, ### A single file node only tracks the file at self.prefix
#     _output=[], 
#     ):
#     '''
#     #### One can also use directly move the output file, but this would break the upstream integrity 
#     #### and is hence not recommended
#     '''
#     output = prefix
#     haso = os.path.exists(output)
#     hasi = os.path.exists(input)

#     shutil.copy2(input, self.prefix+'.temp')
#     shutil.move(self.prefix +'.temp', self.prefix)

def backup(self,prefix):
	key = 'subflow..random_seq..output..seq'
	# self.runner(LinkEvent, prefix)
	self.runner(LinkEvent, prefix+'.' + key, resolve_spiper(flow,key))

def get_fasta(self, prefix,
	_depends = [Depend('bin://curl'),Depend('bin://gzip')],
	_resp = spiper.types.HttpResponseContentHeader('https://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Wuhan_seafood_market_pneumonia_virus/bigZips/chromFa.tar.gz'),
	_output = ['fasta','cmd']):
	with (self.prefix_named/'_temp').makedirs_p() as d:
		CMD = ['curl','-LC0',_resp.url,
		'|','tar','-xvzf-',]
		stdout = spiper.types.LoggedShellCommand(CMD)
		res = d.glob('*.fa')
		assert len(res)==1
		res[0].move(self.output.fasta)
	d.rmtree_p()



def get_genepred(self,prefix,
	_resp = spiper.types.HttpResponseContentHeader('https://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Wuhan_seafood_market_pneumonia_virus/database/ncbiGene.txt.gz'),
	_IMAGE = Depend('docker://quay.io/biocontainers/ucsc-genepredtogtf:377--h35c10e6_2'),
	_output = ['genepred','gtf','cmd'],
	):
	CMD = ['curl','-LC0',_resp.url,
	'|','gzip -d | cut -f2- >',self.output.genepred,
	]

	LoggedShellCommand(CMD,self.output.cmd,mode='w')
	CMD = ['genePredToGtf','file',self.output.genepred, self.output.gtf]
	LoggedSingularityCommand(self.prefix_named, CMD, _IMAGE, self.output.cmd,mode='a')

@Flow
def test_job(self,prefix, _THREADS= 2,
	_dummy_dep = RPO('spiper_mock_flow@https://github.com/shouldsee/spiper_mock_flow/tarball/master','TOPLEVEL:workflow'),
	_output=[]):
	curr = self.runner(get_fasta, prefix)
	curr = self.runner(get_genepred, prefix)
	curr = self.runner(workflow, 
		prefix+'.sample1', 
		prefix+'.hisat2/wuhan-ncov19' , 
		self.subflow['get_fasta'].output.fasta,
		self.subflow['get_genepred'].output.gtf,
		'./test_data/test_R1_.fastq',
		'./test_data/test_R2_.fastq',
		_THREADS
		)
	return self

if __name__ == '__main__':
	from pprint import pprint
	from spiper.runner import force_run,cache_run
	cache_run(test_job,'_temp_build/root')
	'python3 -m spiper run     $PACKAGE TOPLEVEL:test_job _temp_build/root'

######
if 0:
	################################### TBC afterwards ############################



	def get_htseq():
		_ = '''
		#### htseq-count is too slow and not used
		htseq-count 
		-s reverse 
		-f bam 809_S1.bam 
		/home/feng/ref/Arabidopsis_thaliana_TAIR10/annotation/genes.gtf 
		-r pos -o 809_S1.htseq.sam >809_S1.htseq.count 2>809_S1.htseq.log
		'''
		pass

