#!/usr/bin/env python3

import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Pulls all updates from external repositories.')
parser.add_argument('-b', '--branch', help='Git branch of GMCORE to pull', default='master')
args = parser.parse_args()

project_root = os.getcwd()

def run(cmd):
	print(f'==> {cmd}')
	res = subprocess.run(cmd, shell=True, check=True)

def pull(path):
	print(f'[Notice]: Check updates of {path}.')
	os.chdir(path)
	run('git checkout master')
	run('git pull')
	os.chdir(project_root)

pull('lib/container')
pull('lib/container')
pull('lib/datetime')
pull('lib/fiona')
pull('lib/flogger')
pull('lib/string')

print(f'[Notice]: Check updates of GMCORE.')
run(f'git pull origin {args.branch}')
