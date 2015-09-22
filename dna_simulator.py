import random
import numpy as np
import os
import errno
import argparse

amino_acids = ['A','C','G','T']

def gen_sequence(N):
    return ''.join([random.choice(amino_acids) for _ in range(N)])
    
def insert_alu(seq,alu,n):
	seq = list(seq)
	alu = list(alu)
	alu_pos = []
	for i in range(n):
		alu_pos.append(random.randrange(len(seq)+i*n))
	alu_pos.sort()
	for p in alu_pos:
		seq = seq[:p] + alu + seq[p:]
	return ''.join(seq), alu_pos
	
def alter_sequence(seq, freq):
	alt_seq = list(seq)[:]
	for i in range(len(seq)):
		if (random.random() < freq):
			alt_seq[i] = random.choice(amino_acids)
	return ''.join(alt_seq)
	
def gen_pair_reads(seq, L, coverage, error, mean_dist, std_dev):
	N = int((coverage * len(seq))//(L*2))
	seq = list(seq)
	read_pairs = []
	for _ in range(N):
		d = int(round(random.normalvariate(mean_dist,std_dev)))
		l_pos = random.randrange(len(seq)-(2*L+d))
		r_pos = l_pos + L + d
		l_read = seq[l_pos:l_pos+L]
		r_read = seq[r_pos:r_pos+L]
		if(random.random() < error):
			l_read[random.randrange(L)] = random.choice(amino_acids)
		if(random.random() < error):
			r_read[random.randrange(L)] = random.choice(amino_acids)
		read_pairs.append((l_read,r_read))
	return read_pairs
	
def create_test_dir(test_dir):
    try:
        os.makedirs(test_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("N", help="Give the input size for the length of the DNA string to generate.", type=int)
	parser.add_argument("Alu", help="The number of Alu insertions to make in the donor sequence.", type=int)
	parser.add_argument("-c", "--coverage", help="The sampling coverage. Should be between 10 and 50.", type=int, default=30)
	parser.add_argument("-L", "--read_length", help="The length of reads generated in the simulator")
	args = parser.parse_args()
	N = args.N
