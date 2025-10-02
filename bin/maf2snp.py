#!/usr/bin/env python
import argparse
from itertools import zip_longest


def parse_args():
    # Argument definition
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--maf",
        metavar="maffile.maf",
        type=str,
        help="Input alignments",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out",
        metavar="prefix",
        type=str,
        help="Output file prefix.",
        required=False,
        default="snps",
    )
    parser.add_argument(
        "-d",
        "--delim",
        metavar="char",
        type=str,
        help="Delimiter.",
        required=False,
        default=",",
    )
    parser.add_argument(
        "--genomes",
        type=str,
        help="Genomes to consider.",
        nargs="+",
        required=True,
    )
    return parser.parse_args()


def reader(maf):
    lines = []
    infile = open(maf)
    for line in infile:
        if "#" in line:
            continue
        elif line[0] == "a":
            continue
        elif line[0] == "s":
            lines.append(line)
            continue
        elif line.strip() == "" and len(lines) != 0:
            yield lines
            lines = []
        else:
            lines = []
    if len(lines) > 0:
        yield lines


SFX = {
    ",": "csv",
    ";": "csv",
    "\t": "tsv",
    " ": "txt",
}


def get_matching_prefix(string, prefixes):
    for prefix in prefixes:
        if string.startswith(prefix):
            return prefix
    raise Exception(f'Prefix not found for: {string}')


def main():
    # Get arguments
    args = parse_args()
    # Define unique list of input genomes
    genomes = set(args.genomes)

    # Prepare output
    if args.out == '-' or args.out == '/dev/stdout':
        ofname = '/dev/stdout'
    else:
        ofname = f"{args.out}.{SFX.get(args.delim, 'txt')}"
    outfile = open(ofname, "w")

    # Process the data
    outfile.write("CHROM\tPOS\t{}\n".format("\t".join(args.genomes)))
    for n, aln in enumerate(reader(args.maf)):
        _, chrxspp, bpi, alnlen, _, _, _ = aln[0].split()
        try:
            _, chromosome = chrxspp.split(".", 1)
        except:
            raise Exception(f"Impossible to split: {chrxspp}")
        splits = [a.split() for a in aln]            
        sequences = {get_matching_prefix(a[1], args.genomes): a[-1] for a in splits}
        for genome in genomes.difference(set(sequences.keys())):
            sequences[genome] = "-" * int(alnlen)
        zipped = list(zip_longest(*[sequences[g] for g in args.genomes]))
        for n, bp in enumerate(range(int(bpi), int(bpi) + int(alnlen))):
            outfile.write("{}\t{}\t{}\n".format(chromosome, bp, "\t".join(zipped[n])))


if __name__ == "__main__":
    main()
