from config import *
from utils import load_msa_to_list, get_positional_counts, find_node_weights, \
					find_edge_correlations
from analyze_utils import count_num_edges_per_seq, gen_histogram

# load MSA alignment output to 2D list for each seqid/seq
msa_seqids, msa_seqs, msa_num_seqs, msa_seq_len = load_msa_to_list(class_msa_loc)

# generate count dict of kmers at each position
msa_position_counts = get_positional_counts(msa_seqs, characters_of_interest, msa_seq_len)

# We consider a node to be each observed unique nucleotide at each position of the MSA
observed_nodes_w_weights = find_node_weights(msa_position_counts, characters_of_interest, \
							msa_seq_len, msa_num_seqs, \
							class_node_storage)

# We consider each edge to be each observed occurence of a uniq nt at 
# one position to another position in a sequence
# Get p-values of each edge via Fishers Exact Test
observed_edges_w_weights = find_edge_correlations(msa_seqs, characters_of_interest, \
							msa_seq_len, msa_num_seqs, \
							msa_position_counts, \
							class_edge_storage)

# now conduct classification component of analysis
new_seq_tracker, family_seq_tracker = count_num_edges_per_seq(msa_seqids, msa_seqs, msa_seq_len, \
											observed_edges_w_weights, \
											characters_of_interest, class_max_p_value)

# create distribution
gen_histogram(newseq_test_name, histo_out_dir, \
					new_seq_tracker, family_seq_tracker)