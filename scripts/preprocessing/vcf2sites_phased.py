#!/usr/bin/env python3

import gzip
import re
import sys

helpMsg = '''
Usage:
  python3 vcf2sites_phased.py input.vcf.gz output_prefix [chrom start end]

Examples:
  # Process entire VCF
  python3 vcf2sites_phased.py mydata.vcf.gz output_prefix

  # Process specific region (chromosome 2, positions 100000–200000)
  python3 vcf2sites_phased.py mydata.vcf.gz 1 100000 200000
'''

def parseVCFline(vcf_line, gt_map):
    fields = vcf_line.strip().split(b'\t')
    chrom = fields[0].decode()
    pos = int(fields[1])

    bases = []
    for f in fields[9:]:
        gt = f.decode().split(":")[0]  # pega apenas GT
        if '|' not in gt:
            bases.extend(['N', 'N'])  # não faseado? pula
            continue
        a1, a2 = gt.split('|')
        try:
            b1 = gt_map[int(a1)]
            b2 = gt_map[int(a2)]
            bases.extend([b1, b2])
        except (ValueError, IndexError):
            bases.extend(['N', 'N'])

    return chrom, pos, ''.join(bases).replace('*', 'N')

def is_invariant(seq):
    return seq.count(seq[0]) == len(seq)

def main():
    if len(sys.argv) < 3:
        print(helpMsg)
        sys.exit(1)

    vcf_path = sys.argv[1]
    out_prefix = sys.argv[2]

    region_mode = len(sys.argv) == 6
    if region_mode:
        chrom = sys.argv[3]
        win_start = int(sys.argv[4])
        win_end = int(sys.argv[5])
        print(f"Processing region {chrom}:{win_start}-{win_end}...")
    else:
        print("Processing entire VCF...")

    open_func = gzip.open if vcf_path.endswith('.gz') else open
    pos_to_site = {}
    header_lines = []
    inferred_chrom = None

    with open_func(vcf_path, 'rb') as infile:
        for line in infile:
            if line.startswith(b'##'):
                continue
            elif line.startswith(b'#'):
                samples = line.strip().split(b'\t')[9:]
                hap_names = b'\t'.join(
                    sum(([s + b"_1", s + b"_2"] for s in samples), [])
                )
                header_lines.append(f'NAMES\t{hap_names.decode()}')
            else:
                fields = line.strip().split(b'\t')
                chrom_rec = fields[0].decode()
                pos = int(fields[1])
                ref = fields[3].decode()
                alt = fields[4].decode().split(',') if fields[4] != b'.' else []
                gt_map = [ref] + alt

                if not region_mode and inferred_chrom is None:
                    inferred_chrom = chrom_rec

                if region_mode:
                    if chrom_rec != chrom or not (win_start <= pos <= win_end):
                        continue

                _, pos, base_gts = parseVCFline(line, gt_map)
                pos_to_site[pos] = base_gts

    if not region_mode and pos_to_site:
        chrom = inferred_chrom  # define chrom now
        win_start = min(pos_to_site)
        win_end = max(pos_to_site)

    header_lines.insert(1, f'REGION\t{chrom}\t{win_start}\t{win_end}')
    filtered = [(p, s) for p, s in sorted(pos_to_site.items()) if not is_invariant(s)]

    with open(f"{out_prefix}.sites", 'w') as out:
        for h in header_lines:
            out.write(f"{h}\n")
        for pos, gts in filtered:
            out.write(f"{pos}\t{gts}\n")

    print(f" {len(filtered)} variant sites written to {out_prefix}.sites")


if __name__ == "__main__":
    main()


