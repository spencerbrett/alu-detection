import random
import numpy as np
import matplotlib.pyplot as plt
import re
import time


amino_acids = ['A','C','G','T']

def index_sequence(seq, l):
	dna_hash = [[] for _ in xrange(4**l)]
	for i in range(len(seq)-l):
		dna_hash[encode_seq(seq[i:i+l])].append(i)
	return dna_hash

def find_alu_positions(sequence, alu):
	return [m.start() for m in re.finditer(alu, sequence)]

def gen_reads(seq, L, coverage, error):
	N = int((coverage * len(seq))//L)
	reads = []
	for _ in range(N):
		pos = random.randrange(len(seq)-L)
		read = seq[pos:pos+L]
		if (random.random() < error):
			read[random.randrange(L)] = random.choice(amino_acids)
		reads.append(read)
	return reads
		
def encode_seq(seq):
	seq = list(seq)
	n = len(seq)-1
	val = 0
	for i in seq:
		if i == 'A':
			pass
		elif i == 'C':
			val += 1*(4**n)
		elif i == 'G':
			val += 2*(4**n)
		else:
			val += 3*(4**n)
		n -= 1
	return val

def decode_seq(n):
	seq = []
	digs = 'ACGT'
	while n:
		seq.append(digs[n%4])
		n /= 4
	seq.reverse()
	return ''.join(seq)

def compare_seq(s1, s2, err_n):
	if (len(s1) != len(s2)):
		return False
	errors = 0
	for i in range(len(s1)):
		if (s1[i] != s2[i]):
			errors += 1
		if errors > err_n:
			return False
	return errors <= err_n

def map_read(seq,index,read, L, threshold):
	l = L/3
	seq_len = len(seq)
	# Split the read into 3 segments, assuming the length is divisible by 3
	# and l is exactly 1/3 of the read length. l is also the hash size.
	first, second, third = read[0:l], read[l:2*l], read[2*l:3*l]
	# Get the list of candidate positions related to each segment from the index
	f_pos = index[encode_seq(first)]
	s_pos = index[encode_seq(second)]
	t_pos = index[encode_seq(third)]
	f_maps = [pos for pos in f_pos if (compare_seq(read, seq[pos:pos+L], threshold))]
	s_maps = [pos-l for pos in s_pos if (compare_seq(read, seq[pos-l:pos+2*l], threshold))]
	t_maps = [pos-2*l for pos in t_pos if (compare_seq(read, seq[pos-2*l:pos+l], threshold))]
	return list(set(f_maps+s_maps+t_maps))

def alu_read_filter(pair_reads, L, alu, alu_index, threshold):
	candidate_reads = []
	for left, right in pair_reads:
		if map_read(alu, alu_index, left, L, threshold) or map_read(alu, alu_index, right, L, threshold):
			candidate_reads.append((left,right))
	return candidate_reads

def find_mismatch_pairs(seq, index, pair_reads, L, threshold, mean_dist, std_dev):
	mismatch_pairs = []
	for left, right in pair_reads:
		l_maps = map_read(seq,index,left,L,threshold)
		r_maps = map_read(seq,index,right,L,threshold)
		match = False
		for l_map in l_maps:
			for r_map in r_maps:
				target = l_map + L + mean_dist
				# if the right pair maps to target range, we have a matched pair
				if(abs(r_map-target) <= 4*std_dev):
					match = True
		if (l_maps and r_maps) and match == False:
			mismatch_pairs.append((left,right))
	return mismatch_pairs

def find_anchor(pair_read, alu, alu_index, L, threshold):
	left, right = pair_read
	lmap = map_read(alu, alu_index, left, L, threshold)
	rmap = map_read(alu, alu_index, right, L, threshold)
	# both reads map to the alu
	if(lmap and rmap):
		return None, 0
	# left mapped to the alu, so right is the anchor
	elif lmap:
		return "right",lmap[0]
	elif rmap:
		return "left",rmap[0]
	else:
		return None, 0

def calc_alu_position(l_pos, r_pos, anchor,L,mean_dist):
	# left is anchor, so l_pos is position in genome
	# r_pos is position in alu
	if anchor == "left":
		return l_pos + L + mean_dist - r_pos
	elif anchor == "right":
		return r_pos - mean_dist - L - l_pos+300

def cluster_candidates(candidates, cluster_range):
	clusters = []
	while candidates:
		anchor = candidates[0]
		cluster = [item for item in candidates if abs((item - anchor)) <= cluster_range]
		candidates = [item for item in candidates if abs((item - anchor)) > cluster_range]
		clusters.append(cluster)
	for cluster in clusters:
		if len(cluster) < 10:
			clusters.remove(cluster)
	return clusters
	
def round_mean(l):
	return int(round(sum(l)/len(l)))
	
def alu_insertions_baseline(sequence, alu, pair_reads, L, threshold, mean_dist, std_dev):
	ref_alu_positions = find_alu_positions(sequence,alu)
	index = index_sequence(sequence, L/3)
	mismatch_pairs = find_mismatch_pairs(sequence, index, pair_reads, L, threshold, mean_dist, std_dev)
	alu_index = index_sequence(alu, L/3)
	candidate_positions = []
	for left, right in mismatch_pairs:
		anchor, alu_pos = find_anchor((left,right),alu, alu_index, L, threshold)
		l_pos, r_pos = 0,0
		if anchor == None:
			continue
		if anchor == "left":
			l_pos = map_read(sequence,index,left,L,threshold)[0]
			r_pos = alu_pos
		elif anchor == "right":
			l_pos = alu_pos
			r_pos = map_read(sequence,index,right,L,threshold)[0]
		candidate_positions.append(calc_alu_position(l_pos,r_pos,anchor,L,mean_dist))
	candidate_clusters = cluster_candidates(candidate_positions, std_dev*3)
	alt_alu_positions = map(round_mean, candidate_clusters)
	return ref_alu_positions, alt_alu_positions

def alu_insertions_opt(sequence, alu, pair_reads, L, threshold, mean_dist, std_dev):
	ref_alu_positions = find_alu_positions(sequence, alu)
	index = index_sequence(sequence, L/3)
	alu_index = index_sequence(alu, L/3)
	filtered_pairs = alu_read_filter(pair_reads, L, alu, alu_index, threshold)
	mismatch_pairs = find_mismatch_pairs(sequence, index, filtered_pairs, L, threshold, mean_dist, std_dev)
	alu_index = index_sequence(alu, L/3)
	candidate_positions = []
	for left, right in mismatch_pairs:
		anchor, alu_pos = find_anchor((left,right),alu, alu_index, L, threshold)
		l_pos, r_pos = 0,0
		if anchor == None:
			continue
		if anchor == "left":
			l_pos = map_read(sequence,index,left,L,threshold)[0]
			r_pos = alu_pos
		elif anchor == "right":
			l_pos = alu_pos
			r_pos = map_read(sequence,index,right,L,threshold)[0]
		candidate_positions.append(calc_alu_position(l_pos,r_pos,anchor,L,mean_dist))
	candidate_clusters = cluster_candidates(candidate_positions, std_dev*3)
	alt_alu_positions = map(round_mean, candidate_clusters)
	return ref_alu_positions, alt_alu_positions
	
def proc_files():
	alu = ''
	sequence = ''
	pair_reads = []
	with open('ALU_genomeWALU.txt', 'r') as f:
		alu = f.read().rstrip().split(',')[1]
	with open('ref_genomeWALU.txt','r') as f:
		lines = f.read().splitlines()
		sequence = ''.join(lines[2:])
	with open('reads_genomeWALU.txt','r') as f:
		for line in f.read().splitlines()[2:]:
			left, right = line.split(',')
			pair_reads.append((left[:30], right[:30]))
	return sequence, alu, pair_reads
	
def test(seq_len, alt_freq, L, read_err, coverage, ref_alu_n, priv_alu_n, threshold):
	alu = gen_sequence(300)
	ref_seq, ref_alu_pos = insert_alu(gen_sequence(seq_len), alu, ref_alu_n)
	priv_seq, priv_alu_pos = insert_alu(alter_sequence(ref_seq, alt_freq), alu, priv_alu_n)
	pair_reads = gen_pair_reads(priv_seq,L,coverage,read_err, 300, 15)
	index = index_sequence(ref_seq,L/3)
	candidates = [map_read(ref_seq, index, left, L, threshold) for (left,right) in pair_reads] + [map_read(ref_seq, index, right, L, threshold) for (left, right) in pair_reads]
	mappings = []
	for c in candidates:
		for m in c:
			mappings.append(m)
	print "Number of read pairs: "+str(len(pair_reads))
	print "Number of successful read mappings: "+str(len(mappings))
	print "Alu positions in reference: " + str(ref_alu_pos)
	print "Alu positions in private: " + str(priv_alu_pos)
	x = []
	for m in mappings:
		x += range(m,m+L)
	x = np.array(x)
	hist, bins = np.histogram(x, bins=seq_len)
	width = 0.7 * (bins[1]-bins[0])
	center = (bins[:-1] + bins[1:])/2
	plt.bar(center, hist, align='center', width=width)
	plt.show()
	
def test_ref_alu_discovery(seq_len, num_alus):
	alu = gen_sequence(300)
	seq, alu_pos = insert_alu(gen_sequence(seq_len), alu, num_alus)
	alu_pos_exp = find_alu_positions(seq, alu)
	print "Sequence Length: {0}".format(seq_len)
	print "Number of Alu insertions: {0}".format(num_alus)
	print "Number of Alu discovered: {0}".format(len(alu_pos_exp))

def test_alu_read_filter(seq_len,num_alus,coverage, L, threshold):
	alu = gen_sequence(300)
	seq, alu_pos = insert_alu(gen_sequence(seq_len), alu, num_alus)
	pair_reads = gen_pair_reads(seq, L, coverage, 0.01, 100, 5)
	pair_reads_filtered = alu_read_filter(pair_reads, L, alu, threshold)
	print "Sequence Length: {0}".format(seq_len)
	print "Number of Alu insertions: {0}".format(num_alus)
	print "Number of read pairs: {0}".format(len(pair_reads))
	print "Number of read pairs after filter: {0}".format(len(pair_reads_filtered))
	print "Filter percentage: {0}".format((100.0*len(pair_reads_filtered))/len(pair_reads))
	
def test_mismatch_pairs(seq_len,num_alus,coverage,L, threshold):
	alu = gen_sequence(300)
	ref_seq, ref_alu_pos = insert_alu(gen_sequence(seq_len), alu, num_alus)
	index = index_sequence(ref_seq, L/3)
	alt_seq, alt_alu_pos = insert_alu(alter_sequence(ref_seq, 0.01), alu, 1)
	pair_reads = gen_pair_reads(alt_seq, L, coverage, 0.00, 100, 5)
	mismatches = find_mismatch_pairs(ref_seq, index, pair_reads, L, threshold, 100, 5)
	mappings = 0
	for left, right in pair_reads:
		if map_read(ref_seq,index,left,L,threshold):
			mappings += 1
		if map_read(ref_seq,index,right,L,threshold):
			mappings += 1
	print "Sequence Length: {0}".format(seq_len)
	print "Number of Alu insertions: {0}".format(num_alus)
	print "Number of read pairs: {0}".format(len(pair_reads))
	print "Number of non-mapped reads: {0}".format(2*len(pair_reads)-mappings)
	print "Percentage of reads mapped: {0}".format((100.0*mappings)/(2*len(pair_reads)))
	print "Number of mismatch pairs: {0}".format(len(mismatches))

def test_baseline(seq_len, num_alus, coverage,L, threshold):
	alu = gen_sequence(300)
	ref_seq, ref_alu_pos = insert_alu(gen_sequence(seq_len), alu, num_alus)
	index = index_sequence(ref_seq, L/3)
	alt_seq, alt_alu_pos = insert_alu(alter_sequence(ref_seq, 0.05), alu, 2)
	for i in range(len(alt_alu_pos)):
		alt_alu_pos[i] -= i*300
	pair_reads = gen_pair_reads(alt_seq, L, coverage, 0.01, 100, 5)
	exp_ref_alu, exp_alt_alu=alu_insertions_baseline(ref_seq, alu, pair_reads, L, threshold, 100, 5)
	exp_alt_alu.sort()
	print "Sequence Length: {0}".format(seq_len)
	#print "Reference Alu positions: {0}".format(len(ref_alu_pos))
	#print "Experiment Ref positions: {0}".format(len(exp_ref_alu))
	#print "Private Alu positions: {0}".format(len(alt_alu_pos))
	#print "Experiment Pri Alu positions: {0}".format(len(exp_alt_alu))

def find_min_dist(x, l):
	min_dist = abs(x-l[0])
	for i in l[1:]:
		dist = abs(x-i)
		if dist < min_dist:
			min_dist = dist
	return min_dist
	
def determine_accuracy(actual, experimental):
	error = []
	for elem in experimental:
		error.append(find_min_dist(elem, actual))
	return sum(error)/float(len(error))

def test_opt(seq_len, num_alus, coverage,L, threshold):
	std = 2
	alu = gen_sequence(300)
	ref_seq, ref_alu_pos = insert_alu(gen_sequence(seq_len), alu, num_alus)
	index = index_sequence(ref_seq, L/3)
	alt_seq, alt_alu_pos = insert_alu(alter_sequence(ref_seq, 0.05), alu, 2)
	for i in range(len(alt_alu_pos)):
		alt_alu_pos[i] -= i*300
	pair_reads = gen_pair_reads(alt_seq, L, coverage, 0.0, 100, std)
	exp_ref_alu, exp_alt_alu=alu_insertions_opt(ref_seq, alu, pair_reads, L, threshold, 100, std)
	exp_alt_alu.sort()
	print "Sequence Length: {0}".format(seq_len)
	print "Variance: {0}".format(std)
	print "Accuracy: {0}".format(determine_accuracy(alt_alu_pos, exp_alt_alu))
	#print "Reference Alu positions: {0}".format(len(ref_alu_pos))
	#print "Experiment Ref positions: {0}".format(len(exp_ref_alu))
	#print "Private Alu positions: {0}".format(len(alt_alu_pos))
	#print "Experiment Pri Alu positions: {0}".format(len(exp_alt_alu))


def main():
	opt_times6 = []
	for _ in range(5):
		test_opt(100000, 38, 30, 30, 4)
	
	

main()


