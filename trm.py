#!/usr/bin/env python
import os

workpath = '/u/scratch/y/yushendu/RS110416/Zika/'
for i in range(1,434):
  name =  "%03d.m" %i
  print name
  os.system('rm '+workpath+name)

  
