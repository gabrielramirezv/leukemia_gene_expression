#!/bin/bash

# Downloading sequences from NCBI SRA database
nohup prefetch SRR8615998 && fastq-dump --split-files SRR1234567.sra &
