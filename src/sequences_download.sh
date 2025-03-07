#!/bin/bash

# Activar el ambiente donde esta SRA toolkit
conda activate sra-toolkit

# Descargar el archivo .sra
prefetch SRR8615998 && fastq-dump --split-files --skip-technical SRR1234567.sra

