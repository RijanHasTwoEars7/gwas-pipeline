#!/bin/bash
theNum=$1
vecName=${theNum}
stringRun=$2

make_job ()
{
  echo "
  Log        = /home/lconnelly/GWASparallel/logs/${vecName}.log
  Output     =  /home/lconnelly/GWASparallel/logs/${vecName}_tempFiles.out
  Error    = /home/lconnelly/GWASparallel/logs/${vecName}.error
  getenv=true
  request_memory=16GB
  executable                =/home/lconnelly/miniconda3/envs/louisAutocredential/bin/Rscript
  arguments               =/home/lconnelly/GWASparallel/GWASSetaria.R ${stringRun}
  queue" > /home/lconnelly/GWASparallel/logs/${vecName}_condor.sh
  condor_submit /home/lconnelly/GWASparallel/logs/${vecName}_condor.sh

  #get the basename here
}
make_job

