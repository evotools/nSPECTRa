#!/usr/bin/env python
import argparse
import gzip

from tqdm import tqdm


def parse_args():
    # Argument definition
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        metavar="TXT",
        type=str,
        help="Input txt file with the positions (coordinates are zero-based)",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--sites",
        metavar="TXT",
        type=str,
        help="Input list of sites to keep (coordinates are zero-based)",
        required=False,
    )
    parser.add_argument(
        "--sitefile",
        metavar="TXT",
        type=str,
        help="List of sites that have been found.",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--datafile",
        metavar="prefix",
        type=str,
        help="Output data file name.",
        required=False,
        default="data-file.txt",
    )
    parser.add_argument(
        "-c",
        "--configfile",
        metavar="prefix",
        type=str,
        help="Output config file name.",
        required=False,
        default="config-file.txt",
    )
    parser.add_argument(
        "--targets",
        type=str,
        help="Genomes for the target species.",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "--outgroup1",
        type=str,
        help="Outgroup 1 name.",
        required=True,
    )
    parser.add_argument(
        "--outgroup2",
        type=str,
        help="Outgroup 2 name.",
        required=False,
    )
    parser.add_argument(
        "--outgroup3",
        type=str,
        help="Outgroup 3 name.",
        required=False,
    )
    parser.add_argument(
        "--model",
        type=int,
        help="Model for ESFS.",
        default=0,
        choices=[0, 1, 2],
        required=False,
    )
    parser.add_argument(
        "--nrandom",
        type=int,
        help="Random iteration for the ESFS.",
        default=0,
        required=False,
    )
    return parser.parse_args()


def parse_gzipped_line(line):
    return line.decode().strip().split('\t')


def parse_simple_line(line):
    return line.strip().split('\t')


def main():
    # Run main
    args = parse_args()

    # Define input variables
    n_fields = 0
    targets = set(args.targets)
    outgroup1 = args.outgroup1
    outgroup2 = args.outgroup2
    outgroup3 = args.outgroup3
    if outgroup2 is None and outgroup3 is not None:
        raise Exception("Provided outgroup 3 without outgroup 2")
    n_outgroups = 1
    if outgroup2 is not None:
        n_outgroups = 2
    if outgroup3 is not None:
        n_outgroups = 3
    # Variable indexes instantiate
    tgt_idxs = []
    og1_idx = None
    og2_idx = None
    og3_idx = None
    # Define the values to return
    # as counts of [A, C, G, T]
    cnt_idxs = {'a': 0, 'A': 0, 'c': 1, 'C': 1, 'g': 2, 'G': 2, 't': 3, 'T': 3}
    target_cnt = [0, 0, 0, 0]
    outgroup1_cnt = [0, 0, 0, 0]
    outgroup2_cnt = [0, 0, 0, 0]
    outgroup3_cnt = [0, 0, 0, 0]

    # Load sites
    sites = {}
    if args.sites:
        print("Loading target sites...")
        for site in open(args.sites):
            site = site.strip().split()
            sites[f"{site[0]}_{site[1]}"] = True
        print("Found {} sites".format(len(sites)))

    # define the loader
    if args.input.endswith('.gz'):
        buffer = gzip.open
        parser = parse_gzipped_line
    else:
        buffer = open
        parser = parse_simple_line

    # Output file
    file_counter = 1
    row_counter = 0
    outfile = open(f"{args.datafile}.{file_counter}.txt", 'w')
    sitefile = open(f"{args.sitefile}.{file_counter}.txt", 'w')
    # Prepare progress bar for the file
    # Process the input
    print("Processing input file")
    print(f"Saving to file #: {file_counter}")
    pbar = tqdm(total=1000000)
    for n,line in enumerate(buffer(args.input)):
        line = parser(line)
        # First, we work out the header positions
        if n == 0:
            n_fields = len(line)
            # Prepare the index for the targets
            for tgt in targets:
                if tgt not in line:
                    raise Exception(f"Missing key in header: {tgt}")
                tgt_idxs.append( line.index(tgt) )
            # At least one outgroup is needed
            if outgroup1 not in line:
                raise Exception(f"Missing outgroup in header: {outgroup1}")
            og1_idx = line.index(outgroup1)
            # Add at most three outgroups
            if outgroup2:
                if outgroup2 not in line:
                    raise Exception(f"Missing outgroup in header: {outgroup2}")
                og2_idx = line.index(outgroup2)
            if outgroup3:
                if outgroup3 not in line:
                    raise Exception(f"Missing outgroup in header: {outgroup3}")
                og3_idx = line.index(outgroup3)
            continue
        if sites and not sites.get(f"{line[0]}_{line[1]}", False):
                continue
        # If there are already 1M rows, close the file and proceed to the next file
        if row_counter == 1000000:
            pbar.close()
            print(f"Closing file {file_counter} with 1M sites")
            sitefile.close()
            outfile.close()
            row_counter = 0
            file_counter += 1
            print(f"Saving to file {file_counter}")
            pbar = tqdm(total=1000000)
            outfile = open(f"{args.datafile}.{file_counter}.txt", 'w')
            sitefile = open(f"{args.sitefile}.{file_counter}.txt", 'w')
        # Write site id in the file
        sitefile.write("{}\n".format('\t'.join(line[0:2])))
        while len(line) < n_fields:
            line.append('')
        # Then, we process the rest
        for idx in tgt_idxs:
            base = line[idx]
            if base in cnt_idxs.keys():
                target_cnt[cnt_idxs[base]] += 1
        base = line[og1_idx]
        if base in cnt_idxs.keys():
            outgroup1_cnt[cnt_idxs[base]] += 1
        if og2_idx:
            base = line[og2_idx]
            if base in cnt_idxs.keys():
                outgroup2_cnt[cnt_idxs[base]] += 1
        if og3_idx:
            base = line[og3_idx]
            if base in cnt_idxs.keys():
                outgroup3_cnt[cnt_idxs[base]] += 1
        outfile.write('{} {} {} {}\n'.format(
            ','.join(map(str, target_cnt)),
            ','.join(map(str, outgroup1_cnt)),
            ','.join(map(str, outgroup2_cnt)),
            ','.join(map(str, outgroup3_cnt)),
        ))
        target_cnt = [0, 0, 0, 0]
        outgroup1_cnt = [0, 0, 0, 0]
        outgroup2_cnt = [0, 0, 0, 0]
        outgroup3_cnt = [0, 0, 0, 0]
        row_counter += 1
        pbar.update(1)

    # Final closing of the output files & progress bar
    pbar.close()
    sitefile.close()
    outfile.close()

    # Then, save the configuration file
    print("Save configuration file")
    with open(args.configfile, 'w') as outcfg:
        outcfg.write(f"n_outgroup {n_outgroups}\n")
        outcfg.write(f"model {args.model}\n")
        outcfg.write(f"nrandom {args.nrandom}\n")
    print("All done\n\n")

if __name__ == "__main__":
    main()
