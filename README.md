# Read_Mapping

Read_Mapping.py is a Python script developed as a class assignment for the Programming for Bioinformatics course. Its primary function is to map sequencing reads to a reference genome, a fundamental task in bioinformatics for analyzing high-throughput sequencing data.

## Features

Read Mapping: Aligns sequencing reads to a specified reference genome.

Efficiency: Designed to handle large datasets typical in bioinformatics.

Educational Purpose: Serves as a learning tool for understanding the basics of sequence alignment.

## Requirements

Python 3.x: Ensure that Python 3.x is installed on your system.

## Usage

Clone the Repository:
```
git clone https://github.com/BowerH/Read_Mapping.py.git
cd Read_Mapping.py
```

## Prepare Your Data:

Reference Genome: Obtain the reference genome sequence in FASTA format.

Reads: Prepare the sequencing reads in FASTQ format.

## Run the Script:

Execute the script with the required input files:
```
python read_mapping.py -r reference.fasta -q reads.fastq -o output.sam

```
-r or --reference: Path to the reference genome file.

-q or --query: Path to the file containing sequencing reads.

-o or --output: Path to the output file where results will be saved.

## Output:
The script will generate an output file in SAM (Sequence Alignment/Map) format, detailing the alignment of reads to the reference genome.

