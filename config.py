

msa_loc = "resources/queries.aln"
kmer_size = 2
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
