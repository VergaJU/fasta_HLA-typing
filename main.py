#!/home/jacopo/hla-fasta/.venv/bin/python

import argparse
from sanger_typing import read_input, get_exons

# take 2 required arguments from the command line: filder_path and a list of genes to type. The optional argument is the output file name
# argument -h print the help description
parser = argparse.ArgumentParser(description="Type HLA genes from ab1 files")
parser.add_argument("-i", "--input", help="Path to the folder containing the ab1 files", default="./data/")
parser.add_argument("-g","--genes", help="List of genes to type", nargs="+")
parser.add_argument("-o", "--output", help="Output file name", default="./out/")
args = parser.parse_args()

# create a dictionary with the gene as key and the sequence and quality score as value
read_ab1 = read_input.import_ab1(folder_path=args.input)
df = read_ab1.load_ab1_files()

exons = get_exons.get_exons(ab1=df)
alignment = exons.align_local(seqA=df["2F"]["seq"], seqB=df["2R"]["seq"])
alignment