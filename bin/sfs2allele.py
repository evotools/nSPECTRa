#!/usr/bin/env python
import argparse
import sys

import numpy as np


def parse_args():
    # Argument definition
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--sfs",
        metavar="TXT",
        type=str,
        help="SFS probabilities",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="TXT",
        type=str,
        help="Input file for EST-SFS",
        required=True,
    )
    parser.add_argument(
        "-S",
        "--sites",
        metavar="TXT",
        type=str,
        help="Site positions for the SFS input",
        required=False,
    )
    parser.add_argument(
        "-o",
        "--out",
        metavar="prefix",
        type=str,
        help="Output ancestral allele files.",
        required=False,
        default="-",
    )
    return parser.parse_args()


NTS = ['A', 'C', 'G', 'T']


def get_ancestral(maj_p, maj_a, nodes, counts):
    """Define the ancestral allele."""
    # If major allele is more likely than the others, return it 
    if maj_p > 0.5:
        return maj_a
    else:
        alleles = np.array(NTS)
        counts = np.array(map(int, counts))
        min_allele_idx = (counts != np.max(counts)) & (counts != 0)
        min_alleles = set(alleles[min_allele_idx])
        isecs = set(nodes).intersection(min_alleles)
        # If one Nt only is possible, return it
        if len(isecs) == 1:
            return isecs[0]
        # Otherwise, return the oldest node in the tree
        else:
            return nodes[-1]



def main():
    # Get arguments
    args = parse_args()

    # Get major allele for each site
    majors = {}
    with open(args.input) as inputfile:
        for n, line in enumerate(inputfile):
            allele_cnts = [int(i) for i in line.strip().split()[0].split(',')]
            majors[str(n+1)] = {
                'major': NTS[np.argmax(allele_cnts)],
                'vector' : line.strip().split()[0],
                'chrom': None,
                'pos': None
            }

    # Add site mapping
    if args.sites:
        with open(args.sites) as sitesfile:
            for n, line in enumerate(sitesfile):
                chrom, pos = line.strip().split()
                majors[str(n+1)]['chrom'] = chrom
                majors[str(n+1)]['pos'] = pos

    # Define output
    if args.out == '-' or args.out == '/dev/stdout':
        out = sys.stdout
    else:
        out = open(args.out, 'w')

    # Define the ancestral based on the probs
    nts_combinations = [ NT1+NT2 for NT1 in NTS for NT2 in NTS ]
    with open(args.sfs) as sfs:
        for line in sfs:
            if line[0] == '0':
                continue
            site_n, _, maj_prob, probs = line.strip().split(' ', 3)
            maj_prob = float(maj_prob)
            outgr_probs = tuple(map(float, probs.split(' ')))
            max_og_idx = np.argmax(outgr_probs)
            maj_allele = majors[site_n]['major']
            vect = majors[site_n]['vector']
            nodes = nts_combinations[max_og_idx]
            ancestral = get_ancestral(maj_prob, maj_allele, nodes, vect.split(','))
            if args.sites:
                chrom = majors[site_n]['chrom']
                pos = majors[site_n]['pos']
                out.write(f"{site_n}\t{chrom}\t{pos}\t{maj_prob}\t{outgr_probs[max_og_idx]}\t{maj_allele}\t{vect}\t{nodes}\t{ancestral}\n")
            else:
                out.write(f"{site_n}\t{maj_prob}\t{outgr_probs[max_og_idx]}\t{maj_allele}\t{vect}\t{nodes}\t{ancestral}\n")


if __name__ == "__main__":
    main()
