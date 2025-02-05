import argparse
import subprocess
from magnumopus import SAM

parser = argparse.ArgumentParser(description="Map reads and generate consensus sequence")
parser.add_argument("-1", "--read1", required=True, help="Read1 FASTQ file")
parser.add_argument("-2", "--read2", required=True, help="Read2 FASTQ file")
parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file")
parser.add_argument("-s", "--seq_name", help="Sequence name")
args = parser.parse_args()

def run_minimap2(ref_path, read1_path, read2_path, output_sam):
    command = f"minimap2 -ax sr {ref_path} {read1_path} {read2_path} > {output_sam}"
    subprocess.run(command, shell=True, check=True)

output_sam = "mapped_reads.sam"
run_minimap2(args.ref, args.read1, args.read2, output_sam)
sam = SAM.from_sam(output_sam)
if args.seq_name:
    consensus = sam.consensus(args.seq_name)
    print(f">{args.seq_name}_consensus")
    for i in range(0, len(consensus), 60):
        print(consensus[i:i + 60])
else:
    consensus = sam.best_consensus()
    print(">Best Consensus")
    for i in range(0, len(consensus), 60):
        print(consensus[i:i + 60])


