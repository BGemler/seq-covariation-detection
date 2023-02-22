from utils import downselect_seqs_to_observed_features
from scipy.stats import fisher_exact
import csv


def generate_seqstring_combos(downselected_sequences, position_i, num_combos, max_kmer_num_gaps):
	"""
	"""
	uniq_seq_combos = {}
	for seq in downselected_sequences:
		seq_por = seq[position_i:position_i + num_combos]

		if len(seq_por) < 5 and "-" in seq_por:
			continue
		if len(seq_por) >= 5 and seq_por.count("-") > max_kmer_num_gaps:
			continue

		if seq_por not in uniq_seq_combos:
			uniq_seq_combos[seq_por] = 0
		uniq_seq_combos[seq_por] += 1

	uniq_combo_list = []
	for seq_por in uniq_seq_combos:
		count = uniq_seq_combos[seq_por]
		uniq_combo_list.append([seq_por, count])

	return uniq_combo_list


def gen_multi_symbol_edge(multi_symbol_edge_tracking, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, \
														sequences, num_combos):
	"""
	"""
	cont_table = gen_multi_contingency_table(sequences, position_i, seq_por_i, position_j, seq_por_j, num_combos)
	_, pvalue = fisher_exact(cont_table, alternative = 'greater')

	multi_symbol_edge_tracking.append([num_combos, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, pvalue, cont_table])

	return multi_symbol_edge_tracking


def gen_multi_contingency_table(sequences, position_i, seq_por_i, position_j, seq_por_j, num_combos):
	"""
	"""
	num_i_not_j = 0
	num_j_not_i = 0
	num_i_j = 0
	num_not_i_not_j = 0

	for seq in sequences:
		seq_string_i = seq[position_i : position_i + num_combos]
		seq_string_j = seq[position_j : position_j + num_combos]

		if seq_string_i == seq_por_i and seq_string_j == seq_por_j:
			num_i_j += 1
		elif seq_string_i == seq_por_i and seq_string_j != seq_por_j:
			num_i_not_j += 1
		elif seq_string_i != seq_por_i and seq_string_j == seq_por_j:
			num_j_not_i += 1
		elif seq_string_i != seq_por_i and seq_string_j != seq_por_j:
			num_not_i_not_j += 1

	cont_table = [[num_not_i_not_j, num_i_not_j], [num_j_not_i, num_i_j]]

	return cont_table



def multipoint_vectors(msa_seqs, characters_of_interest, \
							msa_position_counts, kmers, \
							msa_num_seqs, msa_seq_len, \
							max_kmer_num_gaps, kmer_strict_p_value, \
							kmer_edge_storage_strict_pvalue, \
							kmer_edge_storage):
	"""
	Calculate correlation values for kmers instead of 1-mers
	"""
	multi_symbol_edge_tracking = []
	for kmer in kmers:
		effective_sequence_length = msa_seq_len - kmer - 1

		# generate sequences indexed by kmer size
		for position_i in range(effective_sequence_length):
			# At each initial position + kmer, find symbols that occur at that position
			expressed_symbols_at_position = []

			for symbol_index_i in range(len(characters_of_interest)):
				count = msa_position_counts[position_i][symbol_index_i]
				if count != 0:
					downselected_sequences = downselect_seqs_to_observed_features(msa_seqs, position_i, characters_of_interest[symbol_index_i])
					uniq_combo_list = generate_seqstring_combos(downselected_sequences, position_i, kmer, max_kmer_num_gaps)
					expressed_symbols_at_position.append(uniq_combo_list)

			for uniq_combo_list in expressed_symbols_at_position:
				for seq_por_i, count_i in uniq_combo_list:
					downselected_sequences = downselect_seqs_to_observed_features(msa_seqs, position_i, characters_of_interest[symbol_index_i])

					# check other positions (if room)
					# preceeding
					if position_i >= kmer:
						for position_j in range(position_i - kmer):
							# counts of seq portions at position j from sequences downselected to position_i's seq portion
							position_j_combo_counts = generate_seqstring_combos(downselected_sequences, position_j, kmer, max_kmer_num_gaps)

							for seq_por_j, count_j in position_j_combo_counts:
								multi_symbol_edge_tracking = gen_multi_symbol_edge(multi_symbol_edge_tracking, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, \
																						msa_seqs, kmer)

					# proceeding
					if effective_sequence_length - position_i - kmer > 0:
						for other_position in range(0, effective_sequence_length - position_i - kmer):
							position_j = position_i + other_position + kmer
							position_j_combo_counts = generate_seqstring_combos(downselected_sequences, position_j, kmer, max_kmer_num_gaps)

							for seq_por_j, count_j in position_j_combo_counts:
								multi_symbol_edge_tracking = gen_multi_symbol_edge(multi_symbol_edge_tracking, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, \
																						msa_seqs, kmer)

	with open(kmer_edge_storage, "w") as f:
		print(kmer_edge_storage)
		out = csv.writer(f)
		out.writerow(["Kmer Size", "Vector 1 Start Position i", "Vector 1 String", \
							"Number of Sequences with Vector 1 at Start Position i", \
							"Vector 2 Start Position j", "Vector 2 String", \
							"Number of Sequences with Vector 2 at Start Position j", \
							"Fisher's Exact P-Value", "Contingency Table"])

		for kmer_size, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, pvalue, cont_table in multi_symbol_edge_tracking:
			position_i = position_i - kmer_size
			position_j = position_j - kmer_size
			out.writerow([kmer_size, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, pvalue, cont_table])
	f.close()

	with open(kmer_edge_storage_strict_pvalue, "w") as f:
		print(kmer_edge_storage_strict_pvalue)
		out = csv.writer(f)
		out.writerow(["Kmer Size", "Vector 1 Start Position i", "Vector 1 String", \
							"Number of Sequences with Vector 1 at Start Position i", \
							"Vector 2 Start Position j", "Vector 2 String", \
							"Number of Sequences with Vector 2 at Start Position j", \
							"Fisher's Exact P-Value", "Contingency Table"])

		for kmer_size, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, pvalue, cont_table in multi_symbol_edge_tracking:
			if pvalue > kmer_strict_p_value:
				continue
			position_i = position_i - kmer_size
			position_j = position_j - kmer_size
			out.writerow([kmer_size, position_i, seq_por_i, count_i, position_j, seq_por_j, count_j, pvalue, cont_table])
	f.close()		

	return