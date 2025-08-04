#!/bin/bash

# atac_preprocessing.sh
# Author: Abdelaziz Awad
# Description: Full ATAC-seq preprocessing workflow
# Usage: bash atac_preprocessing.sh

set -euo pipefail

# -----------------------------
# CONFIGURATION
# -----------------------------
# Assign thread count using command line argument if available, otherwise default to 12
THREADS="${1:-12}"
GENOME_INDEX="mm10"  # Bowtie2 index prefix (already built)
CHROM_SIZES="mm10.chrom.sizes"
PEAK_NAME="KO_vs_WT"
HOMER_GENOME="mm10"
OUTDIR="atac_output"

# Create necessary directories
echo "📁 Creating output directories..."
mkdir -p "$OUTDIR" "$OUTDIR/macs3_output" "$OUTDIR/homer_output" "$OUTDIR/motif_locations_output"

# -----------------------------
# 1. DOWNLOAD & CONVERT SRA
# -----------------------------
echo "📥 Downloading and converting SRA files..."
prefetch SRR13827703 SRR13827711
fasterq-dump SRR13827703 SRR13827711 --split-files --threads "$THREADS"
gzip *.fastq

# -----------------------------
# 2. QC + TRIMMING (fastp)
# -----------------------------
echo "🧹 Running fastp QC and trimming..."
fastp -i SRR13827703_1.fastq.gz -I SRR13827703_2.fastq.gz \
      -o "$OUTDIR/WT_R1.trimmed.fastq.gz" -O "$OUTDIR/WT_R2.trimmed.fastq.gz" \
      --html "$OUTDIR/WT_fastp.html" --thread "$THREADS"
fastp -i SRR13827711_1.fastq.gz -I SRR13827711_2.fastq.gz \
      -o "$OUTDIR/KO_R1.trimmed.fastq.gz" -O "$OUTDIR/KO_R2.trimmed.fastq.gz" \
      --html "$OUTDIR/KO_fastp.html" --thread "$THREADS"

# -----------------------------
# 3. ALIGNMENT (Bowtie2)
# -----------------------------
echo "🧬 Aligning reads with Bowtie2..."
bowtie2 -x "$GENOME_INDEX" -1 "$OUTDIR/WT_R1.trimmed.fastq.gz" -2 "$OUTDIR/WT_R2.trimmed.fastq.gz" -S "$OUTDIR/WT.sam" -p "$THREADS"
bowtie2 -x "$GENOME_INDEX" -1 "$OUTDIR/KO_R1.trimmed.fastq.gz" -2 "$OUTDIR/KO_R2.trimmed.fastq.gz" -S "$OUTDIR/KO.sam" -p "$THREADS"

# -----------------------------
# 4. SORT & INDEX BAM
# -----------------------------
echo "📦 Converting, sorting, and indexing BAM files..."
samtools view -bS "$OUTDIR/WT.sam" | samtools sort -o "$OUTDIR/WT.sorted.bam" --threads "$THREADS"
samtools view -bS "$OUTDIR/KO.sam" | samtools sort -o "$OUTDIR/KO.sorted.bam" --threads "$THREADS"
samtools index "$OUTDIR/WT.sorted.bam"
samtools index "$OUTDIR/KO.sorted.bam"
rm "$OUTDIR/WT.sam" "$OUTDIR/KO.sam" # Clean up intermediate SAM files

# -----------------------------
# 5. FILTER BAM FILES
# -----------------------------
echo "🧯 Filtering for unmapped reads, chrM, and proper pairs..."
samtools view -F 4 -f 2 -q 20 -@ "$THREADS" "$OUTDIR/WT.sorted.bam" | grep -v 'chrM' | samtools view -b > "$OUTDIR/WT.filtered.bam"
samtools view -F 4 -f 2 -q 20 -@ "$THREADS" "$OUTDIR/KO.sorted.bam" | grep -v 'chrM' | samtools view -b > "$OUTDIR/KO.filtered.bam"
samtools index "$OUTDIR/WT.filtered.bam"
samtools index "$OUTDIR/KO.filtered.bam"
rm "$OUTDIR/WT.sorted.bam" "$OUTDIR/KO.sorted.bam" # Clean up intermediate sorted BAMs

# -----------------------------
# 6. REMOVE DUPLICATES (Picard)
# -----------------------------
echo "🧽 Removing duplicates with Picard..."
picard MarkDuplicates I="$OUTDIR/WT.filtered.bam" O="$OUTDIR/WT.cleaned.bam" M="$OUTDIR/WT.metrics.txt" REMOVE_DUPLICATES=true
picard MarkDuplicates I="$OUTDIR/KO.filtered.bam" O="$OUTDIR/KO.cleaned.bam" M="$OUTDIR/KO.metrics.txt" REMOVE_DUPLICATES=true
samtools index "$OUTDIR/WT.cleaned.bam"
samtools index "$OUTDIR/KO.cleaned.bam"
rm "$OUTDIR/WT.filtered.bam" "$OUTDIR/KO.filtered.bam" # Clean up intermediate filtered BAMs

# -----------------------------
# 7. PEAK CALLING (MACS3)
# -----------------------------
echo "🚩 Calling peaks with MACS3..."
macs3 callpeak \
  -t "$OUTDIR/KO.cleaned.bam" \
  -c "$OUTDIR/WT.cleaned.bam" \
  -f BAMPE -g mm \
  -n "$PEAK_NAME" \
  --outdir "$OUTDIR/macs3_output" \
  --nomodel --shift -100 --extsize 200 \
  -B --SPMR --keep-dup all

# -----------------------------
# 8. COVERAGE TRACKS (BigWig)
# -----------------------------
echo "📊 Creating BigWig tracks..."
curl -s -o "$OUTDIR/$CHROM_SIZES" https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
sort -k1,1 -k2,2n "$OUTDIR/macs3_output/${PEAK_NAME}_treat_pileup.bdg" > "$OUTDIR/macs3_output/KO_treat.sorted.bdg"
bedGraphToBigWig "$OUTDIR/macs3_output/KO_treat.sorted.bdg" "$OUTDIR/$CHROM_SIZES" "$OUTDIR/macs3_output/KO_treat.bw"
# Remove intermediate bedGraph file to save space
rm "$OUTDIR/macs3_output/KO_treat.sorted.bdg"

# -----------------------------
# 9. MOTIF ANALYSIS (HOMER)
# -----------------------------
echo "🔎 Running HOMER motif analysis..."
findMotifsGenome.pl "$OUTDIR/macs3_output/${PEAK_NAME}_peaks.narrowPeak" "$HOMER_GENOME" "$OUTDIR/homer_output/" -size given
annotatePeaks.pl "$OUTDIR/macs3_output/${PEAK_NAME}_peaks.narrowPeak" "$HOMER_GENOME" > "$OUTDIR/annotated_peaks.txt"

# Extra motif mapping
findMotifsGenome.pl "$OUTDIR/macs3_output/${PEAK_NAME}_peaks.narrowPeak" "$HOMER_GENOME" "$OUTDIR/motif_locations_output/" -find "$OUTDIR/homer_output/homerMotifs.all.motifs" > "$OUTDIR/motif_locations.txt"

# -----------------------------
# 10. MERGE MOTIF LOCATIONS WITH PEAKS
# -----------------------------
echo "📎 Merging motif info with peaks..."
# Use a more robust and streamlined joining method
join -t $'\t' \
    <(tail -n +2 "$OUTDIR/motif_locations.txt" | sort -k1,1) \
    <(awk 'BEGIN {OFS="\t"} {print $4, $0}' "$OUTDIR/macs3_output/${PEAK_NAME}_peaks.narrowPeak" | sort -k1,1) \
    > "$OUTDIR/merged_motifs_and_peaks.tsv"

# Create a proper header and prepend it
(echo -e "PeakID\tOffset\tSequence\tMotifName\tStrand\tMotifScore\tchr\tstart\tend\tname\tscore\tstrand\tsignal\tpval\tqval\tsummit" && cat "$OUTDIR/merged_motifs_and_peaks.tsv") > "$OUTDIR/final_merged_motifs.tsv"

# Clean up intermediate files
rm "$OUTDIR/merged_motifs_and_peaks.tsv"

# -----------------------------
# DONE
# -----------------------------
echo "✅ ATAC-seq preprocessing complete. All outputs are in the '$OUTDIR' directory."
