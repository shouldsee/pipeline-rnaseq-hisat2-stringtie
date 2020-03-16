#!/usr/bin/env python
#from setuptools import setup
### dummy touch
from distutils.core import setup
import os,glob,sys
assert sys.version_info >= (3,5),('Requires python>=3.5, found python==%s'%('.'.join([str(x) for x in sys.version_info[:3]])))

config = dict(
	name='pipeline_rnaseq_hisat2_stringtie',
	version = '0.0.1',
	 packages=['.'],
)

dict(
	license='MIT',
	author='Feng Geng',
	author_email='shouldsee.gem@gmail.com',
	# long_description=open('README.md').read(),
	classifiers = [
	'Programming Language :: Python :: 3.5',
	'Programming Language :: Python :: 3.7',
	],
	install_requires=[
	],
)
setup(**config)