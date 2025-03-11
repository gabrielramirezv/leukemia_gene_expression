#!/bin/bash

# Descargar las secuencias de RNA-Seq del repositorio SRA en segundo plano
nohup prefetch SRR8615998 &

# Convertir el archivo .sra a formato FASTQ, separando los pares de lecturas y omitiendo datos técnicos
fastq-dump --split-files --skip-technical SRR8615998/SRR8615998.sra

# Evaluar la calidad de las lecturas con FastQC
fastqc SRR8615998_1.fastq
fastqc SRR8615998_2.fastq

# Revisar los reportes HTML generados para verificar la calidad de las secuencias
# Como la calidad de las secuencias es buena, no es necesario aplicar trimming con trim_galore o trimmomatic

# Construir un índice de transcritos para Kallisto a partir de la referencia de Homo sapiens
kallisto index -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx Homo_sapiens.GRCh38.cdna.all.fa

# Cuantificar la expresión génica usando Kallisto con las lecturas obtenidas
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o alineamiento_SRR8615998 SRR8615998_1.fastq SRR8615998_2.fastq
