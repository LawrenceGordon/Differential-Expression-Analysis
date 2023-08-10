#! /usr/bin/env python3

import argparse
import pandas as pd

# parses command-line arguments with paladin tsvs and outname
def parse_args():
    parser = argparse.ArgumentParser(description="Combines paladin TSVs into count matrix")
    parser.add_argument("pal", nargs="+", help="paladin tsvs")
    parser.add_argument("--minqual", "-q", default=0.0, type=float, help="minimum average quality to filter by")
    parser.add_argument("--out", "-o", default="paladin_matrix.txt", help="output name for final file")
    return parser.parse_args()

# parse paladin tsvs into dictionaries
def parse_paladin_tsv(tsv_files, min_qual=0.0):
    pal_tsv_dicts = {}

    for pal_tsv in tsv_files:
        with open(pal_tsv) as pal:
            curr_pal = {}
            for line in pal:
                if line.startswith("Count"):
                    continue
                line = line.rstrip().split("\t")
                count = line[0]; avg_qual = line[2]; uniprot_id = line[4]
                if float(avg_qual) < min_qual:
                    continue
                curr_pal[uniprot_id] = count
        pal_tsv_dicts[pal_tsv] = curr_pal
    
    return pal_tsv_dicts

# combine paladin tsvs into matrix
def combine_tsv_to_matrix(pal_tsv_dicts, out_file):
    samples = sorted(list(pal_tsv_dicts.keys()))

    # get master list of all paladin alignments
    all_uniprot_ids = sorted(list(set().union(*[pal_tsv.keys() for pal_tsv in list(pal_tsv_dicts.values())])))

    with open(out_file, "w") as out:
        out.write("\t{0}\n".format(("\t").join(samples)))
        for id in all_uniprot_ids:
            out_list = [id]
            for sample in samples:
                try:
                    out_list.append(pal_tsv_dicts[sample][id])
                except KeyError:
                    out_list.append("0")
            out.write("{0}\n".format(("\t").join(out_list)))
            
def main():
    args = parse_args()
    pal_tsv_dicts = parse_paladin_tsv(args.pal, args.minqual)
    combine_tsv_to_matrix(pal_tsv_dicts, args.out)
    
main()