import os
import sys

FASTA = os.path.join(os.path.dirname(__file__), "..", "v2", "opt", "human", "hg38", "3utr.fasta")
TARGET = "hg38_knownGene_ENST"

def count_enst(fasta_path):
    count = 0
    with open(fasta_path, "r") as f:
        for line in f:
            if TARGET in line:
                count += 1
    return count

if __name__ == "__main__":
    fasta_path = sys.argv[1] if len(sys.argv) > 1 else FASTA
    count = count_enst(fasta_path)
    print("Count of lines containing '{}': {}".format(TARGET, count))
