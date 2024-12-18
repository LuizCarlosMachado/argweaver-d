# This script converts a haploid VCF file (panTro4.vcf.gz) into a diploid VCF format

import gzip
with gzip.open('panTro4.vcf.gz', 'rt') as inVcf:
    with gzip.open('panTro4_dip.vcf.gz', 'wt') as outVcf:
        outVcf.write('##fileformat=VCFv4.1\n')
        for line in inVcf:
            if line.startswith('#'):
                outVcf.write(line)
            else:
                #outVcf.write(line[:-5]+'1|'+line[-5:])
                outVcf.write(line[:-9]+'\t1|1\n')
