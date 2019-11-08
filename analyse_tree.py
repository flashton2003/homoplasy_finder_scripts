import glob
import dendropy as dpy
from collections import Counter

'''
1. read in tree
2. identify chinese sub-clades
3. count the occurences of homoplasies in each of the chinese sub-clades
4. compare this with the total number of occurences of those homoplasies
'''

def get_node_labels(tree):
	homoplasies = []
	for nd in tree.postorder_node_iter():
		if nd.is_internal():
			nd_annot = nd.annotations.values_as_dict()
			# print nd_annot
			if 'label' in nd_annot:
				homoplasies.extend(nd_annot['label'].split('-'))
	return homoplasies

def print_results(all_homoplasies):
	output_dict = {}
	# print all_homoplasies
	master_list = set([x for sublist in all_homoplasies for x in sublist])
	# print master_list
	homoplasy_counts = []
	for homoplasies in all_homoplasies:
		homoplasy_counts.append(dict(Counter(homoplasies)))
	# print homoplasy_counts
	for h in master_list:
		output_dict[h] = []
		for hcounter in homoplasy_counts:
			if h in hcounter:
				output_dict[h].append(hcounter[h])
			else:
				output_dict[h].append(0)
	for x in output_dict:
		print '\t'.join(map(str, [x, output_dict[x][0], output_dict[x][1], output_dict[x][2]]))





def main(root_dir):
	trees = sorted(glob.glob('%s/*nxs' % root_dir))
	all_homoplasies = []
	for tree_handle in trees:
		# print tree_handle
		tree = dpy.Tree.get_from_path(tree_handle, 'nexus')
		homoplasies = get_node_labels(tree)
		all_homoplasies.append(homoplasies)
		# print
	print_results(all_homoplasies)

	

root_dir = '/Users/flashton/Dropbox/mtb/all_tb_in_asia/homoplasy_finder/results/2019.11.05/post-filipino/sub-trees/'

if __name__ == '__main__':
	main(root_dir)