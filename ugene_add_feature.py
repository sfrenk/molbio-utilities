#!/usr/bin/env python3

# Add feature to Ugene .gb file

import argparse
from Bio import SeqIO
import re
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna

parser = argparse.ArgumentParser(description = "Add feature to Ugene .gb file")
parser.add_argument("file", help = "input .gb file")
parser.add_argument("-g", "--gff", help = "gff3 file containing features to add")

args = parser.parse_args()

infile = SeqIO.read(args.file, "gb")
infile.seq.alphabet= generic_dna

# Get chromosomal coords of Ugene file from the ID
chrom = re.search("([XIV]+)", infile.id).group(1)
start = int(re.search(":([0-9]+),", infile.id).group(1))
end = int(re.search(",([0-9]+)", infile.id).group(1))

# Get features
gff = open(args.gff, "r")
for line in gff:

	feature_chrom = line.strip().split("\t")[0]
	feature_group = line.strip().split("\t")[1]
	feature_start = int(line.strip().split("\t")[3])
	feature_end = int(line.strip().split("\t")[4])
	feature_strand = line.strip().split("\t")[6]
	feature_name = re.search("name=([^;]+);", line.strip().split("\t")[8]).group(1)

	if feature_chrom != chrom or feature_start < start or feature_end > end:
		print("Feature out of range - skipping")
	
	else:
		annot_start = (feature_start-start)
		annot_end = (feature_end-start)

		if feature_strand == "-":
			annot_strand = -1
		else:
			annot_strand = 1

		feature = SeqFeature(FeatureLocation(start = (annot_start), end = (annot_end) + 1), type = 'misc_feature', qualifiers = {"label" : feature_name, "ugene_name" : feature_name, "ugene_group" : feature_group}, strand = annot_strand)
		infile.features.append(feature)

gff.close()

output_file =  open(args.file, "w")
try:
	SeqIO.write(infile, output_file, 'genbank')
except:
	pass
output_file.close()
