#!/bin/bash

# Directorio de trabajo
WORKDIR=$(pwd)

# Crear directorios necesarios
mkdir -p fastq kallisto results fastqc_reports

# Leer muestras desde un archivo
if [ ! -f samples.txt ]; then
    echo "Error: El archivo samples.txt no existe."
    exit 1
fi

# Procesar cada muestra listada en samples.txt
while read -r SAMPLE; do
    echo "Procesando muestra: $SAMPLE"
    
    # Eliminar archivos de bloqueo si existen
    rm -f "$SAMPLE/$SAMPLE.sra.lock"
    
    # Descargar datos
    prefetch "$SAMPLE"
    if [ ! -f "$SAMPLE/$SAMPLE.sra" ]; then
        echo "Error: No se pudo descargar $SAMPLE"
        continue
    fi
    
    # Convertir a formato FASTQ
    fastq-dump --split-files --skip-technical "$SAMPLE/$SAMPLE.sra" -O fastq/
    if [ ! -f "fastq/${SAMPLE}_1.fastq" ] || [ ! -f "fastq/${SAMPLE}_2.fastq" ]; then
        echo "Error: No se generaron archivos FASTQ para $SAMPLE"
        continue
    fi

    # Ejecutar FastQC
    fastqc fastq/${SAMPLE}_1.fastq -o fastqc_reports/
    fastqc fastq/${SAMPLE}_2.fastq -o fastqc_reports/

    # Evaluar calidad con FastQC
    LOW_QUALITY=false
    for REPORT in fastqc_reports/${SAMPLE}_*_fastqc.zip; do
        if [ ! -f "$REPORT" ]; then
            echo "Error: No se generó el reporte FastQC para $SAMPLE"
            continue
        fi
        unzip -q -d fastqc_reports/ "$REPORT"
        SUMMARY_FILE=$(echo "$REPORT" | sed 's/.zip/_fastqc/;s/fastqc_reports\///')/summary.txt
        if [ -f "$SUMMARY_FILE" ] && grep -q "FAIL" "$SUMMARY_FILE"; then
            LOW_QUALITY=true
        fi
    done

    # Si hay problemas de calidad, ejecutar Trimmomatic
    if [ "$LOW_QUALITY" = true ]; then
        echo "Muestra $SAMPLE con baja calidad. Ejecutando Trimmomatic..."
        trimmomatic PE -phred33 \
            fastq/${SAMPLE}_1.fastq fastq/${SAMPLE}_2.fastq \
            fastq/${SAMPLE}_1_trimmed.fastq fastq/${SAMPLE}_1_unpaired.fastq \
            fastq/${SAMPLE}_2_trimmed.fastq fastq/${SAMPLE}_2_unpaired.fastq \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    fi

done < samples.txt

# Descargar el transcriptoma de referencia si no existe
cd kallisto
if [ ! -f Homo_sapiens.GRCh38.cdna.all.fa.gz ]; then
    wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
fi

# Crear el índice de Kallisto si no existe
if [ ! -f Homo_sapiens.GRCh38.cdna.all.idx ]; then
    kallisto index -i Homo_sapiens.GRCh38.cdna.all.idx Homo_sapiens.GRCh38.cdna.all.fa.gz
fi
cd "$WORKDIR"

# Ejecutar Kallisto para cada muestra listada en samples.txt
while read -r SAMPLE; do
    echo "Cuantificando muestra: $SAMPLE"
    if [ ! -f "fastq/${SAMPLE}_1_trimmed.fastq" ] || [ ! -f "fastq/${SAMPLE}_2_trimmed.fastq" ]; then
        echo "Error: No hay archivos FASTQ válidos para $SAMPLE, saltando Kallisto"
        continue
    fi
    kallisto quant -i kallisto/Homo_sapiens.GRCh38.cdna.all.idx -o "kallisto/${SAMPLE}" \
        fastq/${SAMPLE}_1_trimmed.fastq fastq/${SAMPLE}_2_trimmed.fastq
    
    # Renombrar y mover resultado
    mv kallisto/${SAMPLE}/abundance.tsv results/abundance_${SAMPLE}.tsv

done < samples.txt

# Limpieza opcional
rm -r SRR*
rm -r kallisto
