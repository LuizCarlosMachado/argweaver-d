#!/usr/bin/env python3

import gzip
import re
import sys

helpMsg = '''
Usage:
  python3 vcf2sites.py input.vcf.gz output_prefix [chrom start end]

Examples:
  # Process entire VCF
  python3 vcf2sites.py mydata.vcf.gz output_prefix

  # Process specific region (chromosome 2, positions 100000–200000)
  python3 vcf2sites.py mydata.vcf.gz output_prefix 2 100000 200000
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
    if len(sys.argv) < 3:
        print(helpMsg)
        sys.exit(1)

    vcf_path = sys.argv[1]
    out_prefix = sys.argv[2]

    # Check if region was given
    region_mode = len(sys.argv) == 6

    if region_mode:
        chrom = sys.argv[3]
        win_start = int(sys.argv[4])
        win_end = int(sys.argv[5])
        print(f"Processing region {chrom}:{win_start}-{win_end}...")
    else:
        print("Processing entire VCF (no region filtering)...")

    if vcf_path.endswith('.gz'):
        vcf_file = gzip.open(vcf_path, 'rb')
    else:
        vcf_file = open(vcf_path, 'rb')

    pos_to_site = {}
    n_haps = None
    header_written = False

    with vcf_file as infile:
        for line in infile:
            if line.startswith(b'##'):
                continue
            elif line.startswith(b'#'):
                samples = line.strip().split(b'\t')[9:]
                hap_names = b"\t".join(
                    sum(([s + b"_1", s + b"_2"] for s in samples), [])
                )
                n_haps = len(hap_names.split(b'\t'))
                header_lines = [f'NAMES\t{hap_names.decode()}']
            else:
                chrom_rec, pos, gts = parseVCFline(line, len(samples))
                if not region_mode:
                    chrom = chrom_rec  # assign first chromosome seen if not specified
                if region_mode:
                    if chrom_rec != chrom:
                        continue
                    if not (win_start <= pos <= win_end):
                        continue
                pos_to_site[pos] = gts  # keep only last occurrence for each POS

    if not region_mode and pos_to_site:
        win_start = min(pos_to_site.keys())
        win_end = max(pos_to_site.keys())

    # Add REGION line expected by ARGweaver
    header_lines.insert(1, f'REGION\t{chrom}\t{win_start}\t{win_end}')

    # Filter invariant sites
    filtered_sites = [(pos, gts) for pos, gts in sorted(pos_to_site.items()) if not is_invariant(gts, n_haps)]

    with open(f"{out_prefix}.sites", 'w') as outfile:
        for header in header_lines:
            outfile.write(f"{header}\n")
        for pos, gts in filtered_sites:
            outfile.write(f"{pos}\t{gts}\n")

    print(f"✅ {len(filtered_sites)} variant sites saved to {out_prefix}.sites")

if __name__ == "__main__":
    main()

