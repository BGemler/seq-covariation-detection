from config import *
from utils import load_msa_to_list, get_positional_counts, find_node_weights, \
					find_edge_correlations
from plotting import plot_3d_covariation


# load MSA alignment output to 2D list for each seqid/seq
msa_seqids, msa_seqs, msa_num_seqs, msa_seq_len = load_msa_to_list(msa_loc)

# generate count dict of kmers at each position
msa_position_counts = get_positional_counts(msa_seqs, characters_of_interest, msa_seq_len)

# We consider a node to be each observed unique nucleotide at each position of the MSA
observed_nodes_w_weights = find_node_weights(msa_position_counts, characters_of_interest, \
							msa_seq_len, msa_num_seqs, \
							node_storage)

# We consider each edge to be each observed occurence of a uniq nt at 
# one position to another position in a sequence
# Get p-values of each edge via Fishers Exact Test
observed_edges_w_weights = find_edge_correlations(msa_seqs, characters_of_interest, \
							msa_seq_len, msa_num_seqs, \
							msa_position_counts, \
							edge_storage)

# make a 3D plot to view covariation sphere
plot_3d_covariation(characters_of_interest, msa_seq_len, \
							max_p_value, very_sig_p_value, \
							circle_radius, cylinder_height, symbol_points_sizeshrink, \
							vector_color, radius_mult, vertical_extension, \
							observed_edges_w_weights, observed_nodes_w_weights)