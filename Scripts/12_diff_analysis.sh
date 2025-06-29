#!/bin/bash

read -p "Enter reference genome FASTA file: " GENOME_FASTA
read -p "Enter reference annotation file (GFF3 or GTF): " ANNOTATION
read -p "Enter directory containing FASTQ RNA_seq files (.fastq.gz): " FASTQ_DIR
read -p "Enter output directory: " OUTDIR
read -p "Enter number of threads to use [default: 4]: " THREADS
THREADS=${THREADS:-4}

mkdir -p "$OUTDIR"


read -p "Enter number of samples: " N_SAMPLES

declare -a SAMPLES
declare -a CONDITIONS

echo "Enter sample names (without _R1/_R2 suffix) and condition for each:"
for (( i=0; i<$N_SAMPLES; i++ ))
do
  read -p "Sample $((i+1)) name: " SAMPLE
  read -p "Sample $((i+1)) condition (e.g. control): " CONDITION
  SAMPLES+=("$SAMPLE")
  CONDITIONS+=("$CONDITION")
done

# indexing
HISAT2_INDEX="$OUTDIR/hisat2_index"

if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
  echo "Building HISAT2 index..."
  hisat2-build "$GENOME_FASTA" "$HISAT2_INDEX"
else
  echo "HISAT2 index already exists. Skipping build step."
fi

# mapping
echo "Starting alignment with HISAT2..."

for SAMPLE in "${SAMPLES[@]}"
do
  echo "Processing $SAMPLE..."

  R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
  R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: Missing FASTQ files for $SAMPLE. Skipping..."
    continue
  fi

  SAM_OUT="$OUTDIR/${SAMPLE}.sam"
  BAM_OUT="$OUTDIR/${SAMPLE}.sorted.bam"

  hisat2 -x "$HISAT2_INDEX" \
         -1 "$R1" -2 "$R2" \
         -S "$SAM_OUT" \
         --rna-strandness RF \
         -p $THREADS

  echo "Converting SAM to sorted BAM..."
  samtools view -@ $THREADS -bS "$SAM_OUT" | samtools sort -@ $THREADS -o "$BAM_OUT"
  samtools index "$BAM_OUT"
  rm "$SAM_OUT"
done

echo "Alignment completed."

# conversion gff3 to gtf
if [[ "$ANNOTATION" == *.gff3 || "$ANNOTATION" == *.gff ]]; then
  echo "Converting GFF3 to GTF..."
  CONVERTED_GTF="$OUTDIR/converted_annotation.gtf"
  gffread "$ANNOTATION" -T -o "$CONVERTED_GTF"
  ANNOTATION="$CONVERTED_GTF"
fi

# featurecounts quantifying
echo "Running featureCounts..."

BAM_FILES=()
for SAMPLE in "${SAMPLES[@]}"
do
  BAM="${OUTDIR}/${SAMPLE}.sorted.bam"
  if [ -f "$BAM" ]; then
    BAM_FILES+=("$BAM")
  fi
done

featureCounts -p -a "$ANNOTATION" -o "$OUTDIR/counts.txt" -T $THREADS "${BAM_FILES[@]}"

echo "Quantification completed."

# --- 4. METADATA for DESeq2 ---
echo -e "sample\tcondition" > "$OUTDIR/conditions.txt"
for i in "${!SAMPLES[@]}"
do
  echo -e "${SAMPLES[$i]}\t${CONDITIONS[$i]}" >> "$OUTDIR/conditions.txt"
done

echo "All steps completed. Results saved in: $OUTDIR"
