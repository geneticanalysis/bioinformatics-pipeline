# Bioinformatics Pipeline for FASTQ Processing

## 🧬 Overview
This bioinformatics pipeline generates synthetic FASTQ data, aligns sequences to a reference genome using BWA, and detects genetic variants with bcftools. It automates the process of:
- **Generating synthetic FASTQ files** with random DNA sequences.
- **Creating a reference genome (FASTA)**.
- **Analyzing FASTQ files** (read count, length distribution).
- **Aligning sequences to the reference genome** using `BWA`.
- **Generating a sorted BAM file** using `Samtools`.
- **Calling variants to create a VCF file** using `bcftools`.

## 🚀 Features
- Generates synthetic **FASTQ** and **FASTA** files.
- Uses **BWA** to align reads to a reference genome.
- Converts **SAM → BAM → Sorted BAM** using **Samtools**.
- Detects genetic variants and outputs a **VCF** file.

## 📦 Dependencies
Ensure you have the following tools installed:
- **Python 3** (with `biopython`, `pandas`)
- **BWA** (for sequence alignment)
- **Samtools** (for BAM processing)
- **bcftools** (for variant calling)

To install Python dependencies, run:
```bash
pip install biopython pandas
