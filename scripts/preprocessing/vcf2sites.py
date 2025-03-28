#!/usr/bin/env python3

import gzip
import re
import sys

helpMsg = '''
Usage:
python3 vcf2filtered_sites.py input.vcf.gz output_prefix chrom start end [--exclude-last]

Optional:
  --exclude-last   Exclude the last sample (e.g., panTro4) from the NAMES line and the haplotype output.
'''

def parseVCFline(vcf_line, n_samples):
    fields = vcf_line.strip().split(b'\t')
    chrom = fields[0].decode()
    pos = int(fields[1])
    ref = fields[3].decode()
    alt_raw = fields[4].decode()
    alt = alt_raw.split(',') if alt_raw != '.' else []
    gt_map = [ref] + alt

    gts = sum([re.split(r'[/|:]', field.decode())[:2] for field in fields[9:]], [])
    def get_base(x):
        if x in {'.', './.'}:
            return 'N'
        try:
            return gt_map[int(x)]
        except (ValueError, IndexError):
            return 'N'
    base_gts = ''.join(map(get_base, gts)).replace('*', 'N')
    return chrom, pos, base_gts

def is_invariant(seq, n_haps):
    return seq.count(seq[0]) == n_haps

def main():
    if not (6 <= len(sys.argv) <= 7):
        print(helpMsg)
        sys.exit(1)

    vcf_path = sys.argv[1]
    out_prefix = sys.argv[2]
    chrom = sys.argv[3]
    win_start = int(sys.argv[4])
    win_end = int(sys.argv[5])
    exclude_last = (len(sys.argv) == 7 and sys.argv[6] == '--exclude-last')

    print(f"Processing region {chrom}:{win_start}-{win_end}...")
    if exclude_last:
        print("Note: the last sample (e.g., panTro4) will be excluded.")

    if vcf_path.endswith('.gz'):
        vcf_file = gzip.open(vcf_path, 'rb')
    else:
        vcf_file = open(vcf_path, 'rb')

    all_sites = []
    n_haps = None
    header_written = False

    with vcf_file as infile:
        for line in infile:
            if line.startswith(b'##'):
                continue
            elif line.startswith(b'#'):
                samples = line.strip().split(b'\t')[9:]
                if exclude_last:
                    samples = samples[:-1]
                hap_names = b"\t".join(
                    sum(([s + b"_1", s + b"_2"] for s in samples), [])
                )
                n_haps = len(hap_names.split(b'\t'))
                header_lines = [
                    f'NAMES\t{hap_names.decode()}',
                    f'REGION\t{chrom}\t{win_start}\t{win_end}'
                ]
            else:
                chrom_rec, pos, gts = parseVCFline(line, len(samples))
                if chrom_rec != chrom:
                    continue
                if win_start <= pos <= win_end:
                    if exclude_last:
                        gts = gts[:-2]
                    all_sites.append((pos, gts))

    # Filter invariant sites
    filtered_sites = [(pos, gts) for (pos, gts) in all_sites if not is_invariant(gts, n_haps)]

    with open(f"{out_prefix}.sites", 'w') as outfile:
        for header in header_lines:
            outfile.write(f"{header}\n")
        for pos, gts in filtered_sites:
            outfile.write(f"{pos}\t{gts}\n")

    print(f"{len(filtered_sites)} variable sites saved to {out_prefix}.sites")

if __name__ == "__main__":
    main()
