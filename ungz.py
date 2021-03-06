#!/usr/bin/env python
import os
import glob
import sys

infile = sys.argv[1]
os.system('gunzip /u/scratch/y/yushendu/Zika/data/'+infile)
