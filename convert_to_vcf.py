from Bio import SeqIO

class Homoplasy():
    def __init__(self, split_line):
        self.snp_sites_position = int(split_line[0])
        self.consistency_index = float(split_line[1])
        self.counts_A = int(split_line[2].split(':')[0])
        self.counts_C = int(split_line[2].split(':')[1])
        self.counts_G = int(split_line[2].split(':')[2])
        self.counts_T = int(split_line[2].split(':')[3])
        self.genomic_position = int
        self.ref_allele = str
        self.alt_alleles = list
        self.multi_allelic = False


def read_in_ss_vcf(ss_vcf_handle):
    with open(ss_vcf_handle) as fi:
        header = []
        variants = []
        for line in fi:
            if line.startswith('#'):
                header.append(line)
            else:
                variants.append(line.split()[1])
    # print len(header)
    return variants

def read_reference_genome(reference_genome_handle):
    seq = SeqIO.read(reference_genome_handle, 'fasta')
    return seq

def read_homoplasy_output(homoplasy_output_handle):
    homoplasies = []
    with open(homoplasy_output_handle) as fi:
        lines = fi.readlines()
        for line in lines[1:]:
            split_line = line.strip().split()
            h = Homoplasy(split_line)
            homoplasies.append(h)
    return homoplasies

def add_genomic_position(homoplasies, variants):
    for h in homoplasies:
        # print h
        ## add in the genomic position of the variant
        h.genomic_position = int(variants[int(h.snp_sites_position) - 1])
        # print h

def add_ref_allele(homoplasies, ref):
    for h in homoplasies:
        h.ref_allele = ref.seq[h.genomic_position - 1]

def check_if_multiallelic(homoplasies):
    for h in homoplasies:
        if len([x for x in [h.counts_A, h.counts_C, h.counts_G, h.counts_T] if x > 0]) > 2:
            # print h.counts_A, h.counts_C, h.counts_G, h.counts_T
            h.multi_allelic = True

def add_alt_alleles(homoplasies):
    for h in homoplasies:
        alts = []
        if h.counts_A > 0:
            alts.append('A')
        if h.counts_C > 0:
            alts.append('C')
        if h.counts_G > 0:
            alts.append('G')
        if h.counts_T > 0:
            alts.append('T')
        if h.ref_allele in alts:
            alts.remove(h.ref_allele)
        h.alt_alleles = alts
        # print len(h.alt_alleles)

def print_vcf(homoplasies):
    vcf = "##fileformat=VCFv4.1\n##contig=<ID=Chromosome,length=4411532>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t16494-ERR3256127\n"
    for h in homoplasies:
        for a in h.alt_alleles:
            vcf += "Chromosome\t%s\t.\t%s\t%s\t.\t.\tHOMOPLASY_FINDER_INDEX=%s\tGT\t1\n" % (h.genomic_position, h.ref_allele, a, h.snp_sites_position)
    print vcf.rstrip('\n')


def main(ss_vcf_handle, homoplasy_output_handle, reference_genome_handle):
    variants = read_in_ss_vcf(ss_vcf_handle)
    ref = read_reference_genome(reference_genome_handle)
    homoplasies = read_homoplasy_output(homoplasy_output_handle)
    add_genomic_position(homoplasies, variants)
    add_ref_allele(homoplasies, ref)
    # print homoplasies[0].genomic_position, homoplasies[0].ref_allele
    check_if_multiallelic(homoplasies)
    add_alt_alleles(homoplasies)
    # print len(homoplasies)
    print_vcf(homoplasies)


'''
0. get snp-sites vcf for the fasta which went into homoplasy finder
1. get homoplasy finder output
2. convert to vcf - change the positions to original genomic positions
    a. convert all the positions to genomic positions
    b. check which base is the reference in the ref genome
    c. check for multi-allelic snps (how to handle these?)
    d. vcf - chromo, pos, ref, alt, etc.
'''

root_dir = '/Users/flashton/Dropbox/mtb/all_tb_in_asia/homoplasy_finder/results/2019.11.05/post-filipino/'
ss_vcf_handle = '%s/2019.11.04_all_L4_masked.snp-sites.vcf' % root_dir
homoplasy_output_handle = '%s/consistencyIndexReport_05-11-19.txt' % root_dir
reference_genome_handle = '/Users/flashton/Dropbox/mtb/reference_genomes/2019.11.05/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa'

if __name__ == '__main__':
    main(ss_vcf_handle, homoplasy_output_handle, reference_genome_handle)