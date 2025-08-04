#!/bin/bash

# atac_preprocessing.sh
# Author: Abdelaziz Awad
# Description: General bulk ATAC-seq preprocessing pipeline
# Usage: bash atac_preprocessing.sh <TREAT_SRA_ID> <CTRL_SRA_ID> <TREAT_NAME> <CTRL_NAME> [THREADS]

set -euo pipefail

# -----------------------------
# INPUTS & CONFIGURATION
# -----------------------------
TREAT_SRA="$1"         # e.g. SRR12345678
CTRL_SRA="$2"          # e.g. SRR87654321
TREAT_NAME="$3"        # e.g. KO
CTRL_NAME="$4"         # e.g. WT
THREADS="${5:-12}"     # default to 12 threads

GENOME_INDEX="mm10"    # Bowtie2 prefix
HOMER_GENOME="mm10"    # HOMER genome
CHROM_SIZES="mm10.chrom.sizes"
PEAK_NAME="${TREAT_NAME}_vs_${CTRL_NAME}"
OUTDIR="atac_output"

# Make folders
mkdir -p "$OUTDIR" "$OUTDIR/macs3_output" "$OUTDIR/homer_output" "$OUTDIR/motif_locations_output"

# -----------------------------
# 1. DOWNLOAD & FASTQ CONVERT
# -----------------------------
echo "📥 Downloading and converting..."
prefetch "$TREAT_SRA" "$CTRL_SRA"
fasterq-dump "$TREAT_SRA" "$CTRL_SRA" --split-files --threads "$THREADS"
gzip "${TREAT_SRA}"*.fastq "${CTRL_SRA}"*.fastq

# -----------------------------
# 2. QC & TRIMMING (fastp)
# -----------------------------
echo "🧹 Trimming adapters with fastp..."
fastp -i "${TREAT_SRA}_1.fastq.gz" -I "${TREAT_SRA}_2.fastq.gz" \
      -o "$OUTDIR/${TREAT_NAME}_R1.trimmed.fastq.gz" -O "$OUTDIR/${TREAT_NAME}_R2.trimmed.fastq.gz" \
      --html "$OUTDIR/${TREAT_NAME}_fastp.html" --thread "$THREADS"

fastp -i "${CTRL_SRA}_1.fastq.gz" -I "${CTRL_SRA}_2.fastq.gz" \
      -o "$OUTDIR/${CTRL_NAME}_R1.trimmed.fastq.gz" -O "$OUTDIR/${CTRL_NAME}_R2.trimmed.fastq.gz" \
      --html "$OUTDIR/${CTRL_NAME}_fastp.html" --thread "$THREADS"

# -----------------------------
# 3. ALIGNMENT (Bowtie2)
# -----------------------------
echo "🧬 Aligning to genome..."
bowtie2 -x "$GENOME_INDEX" -1 "$OUTDIR/${TREAT_NAME}_R1.trimmed.fastq.gz" -2 "$OUTDIR/${TREAT_NAME}_R2.trimmed.fastq.gz" -S "$OUTDIR/${TREAT_NAME}.sam" -p "$THREADS"
bowtie2 -x "$GENOME_INDEX" -1 "$OUTDIR/${CTRL_NAME}_R1.trimmed.fastq.gz" -2 "$OUTDIR/${CTRL_NAME}_R2.trimmed.fastq.gz" -S "$OUTDIR/${CTRL_NAME}.sam" -p "$THREADS"

# -----------------------------
# 4. SORT & INDEX BAM
# -----------------------------
echo "📦 Sorting and indexing BAMs..."
for SAMPLE in "$TREAT_NAME" "$CTRL_NAME"; do
    samtools view -bS "$OUTDIR/$SAMPLE.sam" | samtools sort -@ "$THREADS" -o "$OUTDIR/$SAMPLE.sorted.bam"
    samtools index "$OUTDIR/$SAMPLE.sorted.bam"
    rm "$OUTDIR/$SAMPLE.sam"
done

# -----------------------------
# 5. FILTER chrM + low quality
# -----------------------------
echo "🧯 Filtering BAMs..."
for SAMPLE in "$TREAT_NAME" "$CTRL_NAME"; do
    samtools view -F 4 -f 2 -q 20 -@ "$THREADS" "$OUTDIR/$SAMPLE.sorted.bam" | grep -v 'chrM' | samtools view -b > "$OUTDIR/$SAMPLE.filtered.bam"
    samtools index "$OUTDIR/$SAMPLE.filtered.bam"
    rm "$OUTDIR/$SAMPLE.sorted.bam"
done

# -----------------------------
# 6. REMOVE DUPLICATES
# -----------------------------
echo "🧽 Removing duplicates with Picard..."
for SAMPLE in "$TREAT_NAME" "$CTRL_NAME"; do
    picard MarkDuplicates I="$OUTDIR/$SAMPLE.filtered.bam" O="$OUTDIR/$SAMPLE.cleaned.bam" M="$OUTDIR/$SAMPLE.metrics.txt" REMOVE_DUPLICATES=true
    samtools index "$OUTDIR/$SAMPLE.cleaned.bam"
    rm "$OUTDIR/$SAMPLE.filtered.bam"
done

# -----------------------------
# 7. PEAK CALLING (MACS3)
# -----------------------------
echo "🚩 Peak calling with MACS3..."
macs3 callpeak \
  -t "$OUTDIR/${TREAT_NAME}.cleaned.bam" \
  -c "$OUTDIR/${CTRL_NAME}.cleaned.bam" \
  -f BAMPE -g mm -n "$PEAK_NAME" \
  --outdir "$OUTDIR/macs3_output" \
  --nomodel --shift -100 --extsize 200 \
  -B --SPMR --keep-dup all

# -----------------------------
# 8. COVERAGE TRACKS (BigWig)
# -----------------------------
echo "📊 Generating BigWig tracks..."
curl -s -o "$OUTDIR/$CHROM_SIZES" "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
# Create BigWig for treatment sample
sort -k1,1 -k2,2n "$OUTDIR/macs3_output/${PEAK_NAME}_treat_pileup.bdg" > "$OUTDIR/macs3_output/${TREAT_NAME}_treat.sorted.bdg"
bedGraphToBigWig "$OUTDIR/macs3_output/${TREAT_NAME}_treat.sorted.bdg" "$OUTDIR/$CHROM_SIZES" "$OUTDIR/macs3_output/${TREAT_NAME}_treat.bw"
rm "$OUTDIR/macs3_output/${TREAT_NAME}_treat.sorted.bdg"
# Create BigWig for control sample
sort -k1,1 -k2,2n "$OUTDIR/macs3_output/${PEAK_NAME}_control_lambda.bdg" > "$OUTDIR/macs3_output/${CTRL_NAME}_control.sorted.bdg"
bedGraphToBigWig "$OUTDIR/macs3_output/${CTRL_NAME}_control.sorted.bdg" "$OUTDIR/$CHROM_SIZES" "$OUTDIR/macs3_output/${CTRL_NAME}_control.bw"
rm "$OUTDIR/macs3_output/${CTRL_NAME}_control.sorted.bdg"

# -----------------------------
# 9. MOTIF DISCOVERY (HOMER)
# -----------------------------
echo "🔎 Motif discovery..."
findMotifsGenome.pl "$OUTDIR/macs3_output/${PEAK_NAME}_peaks.narrowPeak" "$HOMER_GENOME" "$OUTDIR/homer_output" -size given
annotatePeaks.pl "$OUTDIR/macs3_output/${PEAK_NAME}_peaks.narrowPeak" "$HOMER_GENOME" > "$OUTDIR/annotated_peaks.txt"

# Optional motif mapping
findMotifsGenome.pl "$OUTDIR/macs3_output/${PEAK_NAME}_peaks.narrowPeak" "$HOMER_GENOME" "$OUTDIR/motif_locations_output" -find "$OUTDIR/homer_output/homerMotifs.all.motifs" > "$OUTDIR/motif_locations.txt"

# -----------------------------
# 10. MERGE MOTIFS + PEAKS
# -----------------------------
echo "📎 Merging motifs with peaks..."
# Use a proper header and join in one step
(echo -e "PeakID\tOffset\tSequence\tMotifName\tStrand\tMotifScore\tchr\tstart\tend\tname\tscore\tstrand\tsignal\tpval\tqval\tsummit"
join -t $'\t' \
    <(tail -n +2 "$OUTDIR/motif_locations.txt" | sort -k1,1) \
    <(awk 'BEGIN{OFS="\t"}{print $4,$0}' "$OUTDIR/macs3_output/${PEAK_NAME}_peaks.narrowPeak" | sort -k1,1)) \
    > "$OUTDIR/final_merged_motifs.tsv"

# -----------------------------
# DONE
# -----------------------------
echo "✅ General ATAC-seq pipeline complete. All results saved to $OUTDIR"
