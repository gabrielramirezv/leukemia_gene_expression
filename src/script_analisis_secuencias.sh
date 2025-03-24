#!/bin/bash
#$ -N TCGA
#$ -cwd
#$ -o salida.out
#$ -e salida.err

# Habilitar opciones estrictas de Bash
set -euo pipefail

# Cargar entorno conda
source /export/apps/bioconda/etc/profile.d/conda.sh

# Crear directorios necesarios si no existen
mkdir -p fastq kallisto results fastqc_reports

# Descargar el transcriptoma de referencia si no existe
cd kallisto
if [ ! -f Homo_sapiens.GRCh38.cdna.all.fa.gz ]; then
    wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
fi

# Descomprimir el archivo si aún está comprimido
if [ -f Homo_sapiens.GRCh38.cdna.all.fa.gz ]; then
    gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
fi

# Crear el índice de Kallisto si no existe
if [ ! -f Homo_sapiens.GRCh38.cdna.all.idx ]; then
    conda run -n kallisto kallisto index -i Homo_sapiens.GRCh38.cdna.all.idx Homo_sapiens.GRCh38.cdna.all.fa
fi
cd -  # Volver al directorio de trabajo

# Función para procesar cada muestra
process_sample() {
    SAMPLE=$(echo "$1" | xargs)  # Eliminar espacios extra
    echo "[$(date +%H:%M:%S)] Procesando muestra: $SAMPLE"

    # Verificar si la muestra ya ha sido procesada
    if [ -f "results/abundance_${SAMPLE}.tsv" ]; then
        echo "   Ya procesada. Saltando..."
        return
    fi

    # Descargar la muestra SRA
    rm -f "$SAMPLE/$SAMPLE.sra.lock"
    conda run -n sra-tools prefetch "$SAMPLE"

    if [ ! -f "$SAMPLE/$SAMPLE.sra" ]; then
        echo "   Error: No se pudo descargar $SAMPLE"
        return
    fi

    # Convertir a formato FASTQ
    conda run -n sra-tools fastq-dump --split-files --skip-technical "$SAMPLE/$SAMPLE.sra" -O fastq/

    if [ ! -s "fastq/${SAMPLE}_1.fastq" ] || [ ! -s "fastq/${SAMPLE}_2.fastq" ]; then
        echo "Error: Archivos FASTQ vacíos o no generados para $SAMPLE"
        return
    fi

    R1="fastq/${SAMPLE}_1.fastq"
    R2="fastq/${SAMPLE}_2.fastq"

    # Ejecutar FastQC para evaluar la calidad
    fastqc "$R1" -o fastqc_reports/
    fastqc "$R2" -o fastqc_reports/

    # Evaluar calidad con FastQC
    LOW_QUALITY=false
    for REPORT in fastqc_reports/${SAMPLE}_*_fastqc.zip; do
        if [ ! -f "$REPORT" ]; then
            echo "Error: No se generó el reporte FastQC para $SAMPLE"
            continue
        fi

        # Extraer y verificar el reporte de calidad
        unzip -q -d fastqc_reports/ "$REPORT"
        SUMMARY_FILE=$(echo "$REPORT" | sed 's/.zip/_fastqc/;s/fastqc_reports\///')/summary.txt

        if [ -f "$SUMMARY_FILE" ] && grep -q "FAIL" "$SUMMARY_FILE"; then
            LOW_QUALITY=true
        fi
    done

    # Si la calidad es baja, ejecutar Trimmomatic para limpiar las secuencias
    if [ "$LOW_QUALITY" = true ]; then
        echo "   Calidad baja. Ejecutando Trimmomatic..."
        trimmomatic PE -threads 4 -phred33 \
            "$R1" "$R2" \
            "fastq/${SAMPLE}_1_trimmed.fastq" "fastq/${SAMPLE}_1_unpaired.fastq" \
            "fastq/${SAMPLE}_2_trimmed.fastq" "fastq/${SAMPLE}_2_unpaired.fastq" \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
        R1="fastq/${SAMPLE}_1_trimmed.fastq"
        R2="fastq/${SAMPLE}_2_trimmed.fastq"
    fi

    # Ejecutar Kallisto para cuantificar la expresión génica
    echo "Cuantificando muestra: $SAMPLE"
    conda run -n kallisto kallisto quant -i kallisto/Homo_sapiens.GRCh38.cdna.all.idx --threads 4 -o "kallisto/${SAMPLE}" "$R1" "$R2"

    # Verificar si se generó el archivo de abundancia
    if [ -f "kallisto/${SAMPLE}/abundance.tsv" ]; then
        mv "kallisto/${SAMPLE}/abundance.tsv" "results/abundance_${SAMPLE}.tsv"
    else
        echo "Error: No se encontró el archivo de abundancia para $SAMPLE"
    fi

    # Limpieza: eliminar archivos temporales específicos de la muestra
    rm -r "$SAMPLE/"
    rm -f "fastq/${SAMPLE}_1.fastq" "fastq/${SAMPLE}_2.fastq"
    rm -f "fastq/${SAMPLE}_1_trimmed.fastq" "fastq/${SAMPLE}_2_trimmed.fastq"
    rm -f "fastq/${SAMPLE}_1_unpaired.fastq" "fastq/${SAMPLE}_2_unpaired.fastq"

    rm -r "fastqc_reports/${SAMPLE}_1_fastqc"
    rm -r "fastqc_reports/${SAMPLE}_2_fastqc"
    rm -r "fastqc_reports/${SAMPLE}_1_*"
    rm -r "fastqc_reports/${SAMPLE}_2_*"

    rm -r "kallisto/${SAMPLE}"

    echo "[$(date +%H:%M:%S)] Finalizado $SAMPLE"
}

# Exportar la función para que pueda ser utilizada por parallel
export -f process_sample

# Verificar que el archivo samples.txt existe
if [ ! -f samples.txt ]; then
    echo "Error: No se encontró samples.txt"
    exit 1
fi

# Leer el archivo samples.txt y procesar las muestras en paralelo
cat samples.txt | grep -v '^$' | parallel -j 12  process_sample