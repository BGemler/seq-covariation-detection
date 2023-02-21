import csv
from scipy.stats import fisher_exact


def load_msa_to_list(msa_loc):
	"""
	This function loads an MSA text output into a 2D list
	where each row of the list is a sequence, gaps included.
	and each column of the list is a position in the MSA
	"""
	seqid_seq_dict = {}
	with open(msa_loc, "r") as f:
		for row in f:
			if row.startswith("CLUSTAL"):
				continue
			if row.replace("\n","").replace(" ","") == "":
				continue
			seqid = row.split()[0]
			row_seq = row.split(seqid)[1].replace("\n","").strip().upper()

			if seqid not in seqid_seq_dict:
				seqid_seq_dict[seqid] = ""
			seqid_seq_dict[seqid] = seqid_seq_dict[seqid] + row_seq
	f.close()

	msa_seqids, msa_seqs = [], []
	for seqid in seqid_seq_dict:
		seq = seqid_seq_dict[seqid]

		msa_seqids.append(seqid)
		msa_seqs.append(seq)

	msa_num_seqs = len(msa_seqids)
	msa_seq_len = len(msa_seqs[0])

	return msa_seqids, msa_seqs, msa_num_seqs, msa_seq_len


def get_positional_counts(msa_seqs, characters_of_interest, msa_seq_len):
	"""
	for each position, computes the counts of each character of interest
	"""
	msa_position_counts = []
	for i in range(msa_seq_len):
		position_counts = []
		for nt in characters_of_interest:
			position_counts.append(0)

		for seq in msa_seqs:
			nt = seq[i]

			if nt not in characters_of_interest:
				continue

			index = characters_of_interest.index(nt)
			position_counts[index] += 1

		msa_position_counts.append(position_counts)

	return msa_position_counts


def find_node_weights(msa_position_counts, characters_of_interest, \
							msa_seq_len, msa_num_seqs, \
							node_storage):
	"""
	obtains a list of every unique combination of nucleotide & position i
	writes out node list for QA/QC
	"""
	observed_nodes_w_weights = []
	for i in range(msa_seq_len):
		position_counts = msa_position_counts[i]
		for j in range(len(characters_of_interest)):
			nt_count_at_pos = position_counts[j]

			if nt_count_at_pos == 0:
				continue

			weight = float(nt_count_at_pos) / float(msa_num_seqs)

			observed_nodes_w_weights.append([i, j, weight])

	with open(node_storage, "w") as f:
		out = csv.writer(f)
		out.writerow(["Position in MSA", "Character of Interest", "Node Weight"])

		for i, j, weight in observed_nodes_w_weights:
			char = characters_of_interest[j]

			out.writerow([i + 1, char, weight])
	f.close()

	return observed_nodes_w_weights


def downselect_seqs_to_observed_features(msa_seqs, i, nt_of_interest):
	"""
	"""
	downselected_seqs = []

	for seq in msa_seqs:
		seq_nt = seq[i]

		if seq_nt == nt_of_interest:
			downselected_seqs.append(seq)

	return downselected_seqs


def generate_contingency_table(msa_seqs, position_i, symbol_i, position_j, symbol_j):
	"""
	"""
	num_i_not_j = 0
	num_j_not_i = 0
	num_i_j = 0
	num_not_i_not_j = 0

	for seq in msa_seqs:
		seq_symbol_i = seq[position_i]
		seq_symbol_j = seq[position_j]

		if seq_symbol_i == symbol_i and seq_symbol_j == symbol_j:
			num_i_j += 1
		elif seq_symbol_i == symbol_i and seq_symbol_j != symbol_j:
			num_i_not_j += 1
		elif seq_symbol_i != symbol_i and seq_symbol_j == symbol_j:
			num_j_not_i += 1
		elif seq_symbol_i != symbol_i and seq_symbol_j != symbol_j:
			num_not_i_not_j += 1

	cont_table = [[num_not_i_not_j, num_i_not_j], [num_j_not_i, num_i_j]]

	return cont_table


def generate_symbol_edge_entry(observed_edges_w_weights, symbol_distribution, characters_of_interest, \
									msa_seqs, pos_1, pos_1_j, pos_2):
	"""
	Generate a confusion matrix for the edge
	Then, compute the Fishers Exact Test to get a p-value
	"""
	pos_1_nt = characters_of_interest[pos_1_j]

	for pos_2_j in range(len(symbol_distribution)):
		if symbol_distribution[pos_2_j] == 0:
			continue

		pos_2_nt = characters_of_interest[pos_2_j]

		cont_table = generate_contingency_table(msa_seqs, pos_1, pos_1_nt, pos_2, pos_2_nt)
		_, pvalue = fisher_exact(cont_table, alternative = 'greater')

		observed_edges_w_weights.append([pos_1, pos_1_j, pos_2, pos_2_j, pvalue, cont_table])

	return observed_edges_w_weights


def find_edge_correlations(msa_seqs, characters_of_interest, \
							msa_seq_len, msa_num_seqs, \
							msa_position_counts, \
							edge_storage):
	"""
	For each edge (calculated deterministically), create confusion matrix
	and calculate Fishers Exact Test pvalue
	"""
	observed_edges_w_weights = []

	# parse over each position of the MSA
	for i in range(msa_seq_len):
		# at position i, find all nts present
		observed_char_indexes_at_i = []
		position_counts = msa_position_counts[i]
		for j in range(len(characters_of_interest)):
			nt_count_at_pos = position_counts[j]

			if nt_count_at_pos > 0:
				observed_char_indexes_at_i.append([j, nt_count_at_pos])

		for j, nt_count_at_pos in observed_char_indexes_at_i:
			# downselect sequnece space to those with nucleotide j at position i
			downselected_seqs = downselect_seqs_to_observed_features(msa_seqs, i, characters_of_interest[j])
			# generate position-count matrix from this downselected list
			downselected_pos_nt_counts = get_positional_counts(downselected_seqs, characters_of_interest, msa_seq_len)
			# Note = we don't use the actual counts here. Just want to know what nucleotides are actually observed 

			# parse over proceeding positions, add edges
			for other_pos in range(1, msa_seq_len - i):
				other_i = other_pos + i
				symbol_distribution = downselected_pos_nt_counts[other_i]
				observed_edges_w_weights = generate_symbol_edge_entry(observed_edges_w_weights, symbol_distribution, characters_of_interest, \
																			msa_seqs, i, j, other_i)

			# parse over preceeding positions, add edges
			for other_pos in range(i):
				other_i = other_pos
				symbol_distribution = downselected_pos_nt_counts[other_i]
				observed_edges_w_weights = generate_symbol_edge_entry(observed_edges_w_weights, symbol_distribution, characters_of_interest, \
																			msa_seqs, i, j, other_i)

	# write out edges for QA/QC
	with open(edge_storage, "w") as f:
		out = csv.writer(f)
		out.writerow(["Position i", "Position i NT", \
						"Position j", "Position j NT", \
						"Fishers Exact PValue", \
						"Confusion Matrix"])

		for pos_1, pos_1_j, pos_2, pos_2_j, pvalue, cont_table in observed_edges_w_weights:
			pos_1_nt = characters_of_interest[pos_1_j]
			pos_2_nt = characters_of_interest[pos_2_j]

			out.writerow([pos_1 + 1, pos_1_nt, pos_2 + 1, pos_2_nt, pvalue, cont_table])
	f.close()

	return observed_edges_w_weights