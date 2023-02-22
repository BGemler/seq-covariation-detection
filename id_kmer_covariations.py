from config import *
from utils import load_msa_to_list, get_positional_counts
from kmer_utils import multipoint_vectors

# load MSA alignment output to 2D list for each seqid/seq
msa_seqids, msa_seqs, msa_num_seqs, msa_seq_len = load_msa_to_list(msa_loc)

# generate count dict of kmers at each position
msa_position_counts = get_positional_counts(msa_seqs, characters_of_interest, msa_seq_len)

multipoint_vectors(msa_seqs, characters_of_interest, \
							msa_position_counts, kmers, \
							msa_num_seqs, msa_seq_len, \
							max_kmer_num_gaps, kmer_strict_p_value, \
							kmer_edge_storage_strict_pvalue, \
							kmer_edge_storage)