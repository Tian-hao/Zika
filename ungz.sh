#!/bin/bash
#call mutations
echo $SGE_TASK_ID
cd /u/scratch/y/yushendu/Zika/split/
gzip `sed -n ${SGE_TASK_ID}p /u/home/y/yushendu/Zika/ref/zipfile.txt`
