from spiper.types import InputFile,OutputFile,File,TempFile, Prefix, Default
from spiper.types import job_result
from spiper.types import Depend
from spiper.runner import list_flatten_strict, job_from_func

from spiper.shell import SingularityShellCommand
from path import Path
import spiper
assert spiper.VERSION >= '0.0.5',spiper.VERSION

def job_trimmomatic(
	self, prefix,
	FASTQ_FILE_1 = InputFile, 
	FASTQ_FILE_2 = InputFile, 
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
		'trimmomatic','PE',
		'-threads', str(THREADS_), 
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
		res = SingularityShellCommand(CMD, _IMAGE, self.output['cmd'])
		return self
		# return job_result( None, CMD, self.output)

def job_hisat2_index(
	self,prefix, 
	FASTA_FILE = InputFile,
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
	res = SingularityShellCommand(CMD, _IMAGE, self.output.cmd)
	return self



def job_hisat2_align(
	self,prefix,
	INDEX_PREFIX = Prefix,
	FASTQ_FILE_1 = InputFile,
	FASTQ_FILE_2 = InputFile,
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
	CMD = [
	 'hisat2',
	 # hisat2_args,
	 '-x', Prefix(INDEX_PREFIX),
	 '-1', File( FASTQ_FILE_1),
	 '-2', File( FASTQ_FILE_2),
	 # '-U', InputFile( FASTQ_FILE_1),
	 # ['-2',InputFile( FASTQ_FILE_2) ] if FASTQ_FILE_2 else [],
	 '-S', File( self.output.bam +'.sam' ),
	 '--threads', str( THREADS_ ),
	 hisat2_args or [
	 '--no-mixed',
	 '--rna-strandness','RF',
	 '--dta',
	 '--fr'],
	 '&>', File( self.output.log),
	]
	res = SingularityShellCommand(CMD, _IMAGE, self.output.cmd)
	# results.append(job_result( None, CMD, self.output))

	_ = '''
	samtools view /home/feng/temp/187R/187R-S1-2018_06_27_14:02:08/809_S1.sam -b --threads 4 -o 809_S1.bam
	'''
	CMD = [	
	'samtools','view',
	File( self.output.bam+'.sam'),
	'--threads',str(THREADS_),
	'-o', 
	File( self.output.bam+'.unsorted'),
	]
	res = SingularityShellCommand(CMD, _IMAGE_SAMTOOLS, self.output.cmd)


	CMD = [
	'samtools','sort',
	File( self.output.bam + '.unsorted'),
	'--threads', str(THREADS_),
	'-o', 
	File( self.output.bam),
	]
	res = SingularityShellCommand(CMD, _IMAGE_SAMTOOLS, self.output.cmd)
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
	res = SingularityShellCommand(CMD, _IMAGE, self.output.cmd)

from spiper.types import Concat
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
		'REMOVE_DUPLICATES=true',
		]
		res = SingularityShellCommand(CMD, _IMAGE, self.output.cmd_log)
		res = SingularityShellCommand(['samtools','index',self.output.bam],
			_IMAGE_SAMTOOLS, 
			self.output.cmd_log,mode='a',
			extra_files = [self.output.bam+'.bai'])



from spiper.types import Flow
@Flow
def workflow(self, prefix, 

	hisat2_cache_prefix = File,
	genome_fasta = File, 
	genome_gtf_file = File,

	fastq1 = File,
	fastq2 = File,

	THREADS_ = int,
	_output=[]
	):
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
	self.runner(
		job_stringtie_count,
		prefix,
		self.subflow['job_picard_dedup'].output.bam,
		genome_gtf_file,
		THREADS_,
		)
	return self

from spiper.types import resolve_spiper
def backup(self,prefix):
	key = 'subflow..random_seq..output..seq'
	self.runner(copy_file, prefix+'.' + key, resolve_spiper(flow,key))

def get_fasta(self, prefix,
	_depends = [Depend('curl'),Depend('gzip')],
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


from spiper.types import LoggedShellCommand
LoggedSingularityCommand = SingularityShellCommand
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
	LoggedSingularityCommand(CMD, _IMAGE, self.output.cmd,mode='a')

@Flow
def test_job(self,prefix, _THREADS= 2,_output=[]):
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

