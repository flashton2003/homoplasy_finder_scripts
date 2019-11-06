from Bio import SeqIO

'''

1. read in alignment
2. read in consistency index findings
3. read in tree
4. select 10 random inconsistent positions
5. take that position from alignment and make a tree annotation file which has that as the tip label
6. look at tree and check that position is inconsistent.
'''

def read_aln(alignment_handle):
	seqs = {}
	seq_iter = SeqIO.parse(alignment_handle, 'fasta')
	for seq in seq_iter:
		# print seq.id
		seqs[seq.id] = str(seq.seq)
	return seqs

def read_in_consistency(consistency_handle):
	inconsistencies = []
	with open(consistency_handle) as fi:
		lines = fi.readlines()
		for line in lines[1:]:
			split_line = line.strip().split()
			inconsistencies.append(int(split_line[0]))
	return inconsistencies

def print_annot(seqs, inconsistencies, output_dir):
	for pos in inconsistencies:
		with open('%s/%s.txt' % (output_dir, pos), 'w') as fo:
			fo.write('\t'.join(['ID', 'pos', '\n']))
			for seq in seqs:
				fo.write('\t'.join([seq, seqs[seq][pos - 1], '\n']))

def main(alignment_handle, consistency_handle, output_dir):
	seqs = read_aln(alignment_handle)
	inconsistencies = read_in_consistency(consistency_handle)
	print_annot(seqs, inconsistencies, output_dir)


root_dir = '/Users/flashton/Dropbox/mtb/all_tb_in_asia/homoplasy_finder/results/2019.11.05/post-filipino'
output_dir = '/Users/flashton/Dropbox/mtb/all_tb_in_asia/homoplasy_finder/results/2019.11.05/post-filipino/annots'
alignment_handle = '%s/2019.11.04_all_L4_masked.snp-sites.fa' % root_dir
consistency_handle = '%s/consistencyIndexReport_05-11-19.txt' % root_dir
# tree_handle = '%s/2019.11.04_all_L4_masked.snp-sites.fa.tree' % root_dir

if __name__ == '__main__':
	main(alignment_handle, consistency_handle, output_dir)