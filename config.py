

msa_loc = "resources/queries.aln"
characters_of_interest = ["A", "C", "T", "G"]

node_storage = "results/observed_nodes.csv"
edge_storage = "results/observed_edges.csv"

# visualization parameters
max_p_value = 0.001
very_sig_p_value = 0.1 * max_p_value

circle_radius = 900.0
cylinder_height = 300.0
vertical_extension = 1.2 # factor for extending the z-axes above the cylinder

vector_color = "white"
symbol_points_sizeshrink = 1 # >1 = smaller circles, <1 = bigger circles for symbols
radius_mult = 5 # scaling for scatter point size


# number of kmers to look for variation in
kmers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
max_kmer_num_gaps = 1
kmer_edge_storage = "kmer_results/observed_kmer_edges.csv"
kmer_strict_p_value = 0.01
kmer_edge_storage_strict_pvalue = "kmer_results/observed_kmer_edges_strict.csv"


# sequence classification
class_msa_loc = "classification_tool/queries_w_new.aln"
class_node_storage = "classification_tool/observed_nodes.csv"
class_edge_storage = "classification_tool/observed_edges.csv"
class_max_p_value = 0.01

histo_out_dir = "classification_tool/histos/"
newseq_test_name = "seq4"