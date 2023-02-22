import matplotlib.pyplot as plt


def gen_histogram(newseq_test_name, histo_out_dir, \
					new_seq_tracker, family_seq_tracker):
	"""
	"""
	plt.hist(family_seq_tracker, bins = 100, label = "Family Sequences", color = 'b')
	plt.axvline(x = new_seq_tracker, color = 'r', label = "New/Test Sequence")
	plt.xlabel("# of Significant Edges / # of Sequence Nodes")
	plt.ylabel("Number of Family Seqs per Bin")
	plt.legend()

	num_family_new_atleast = 0
	for signal in family_seq_tracker:
		if new_seq_tracker >= signal:
			num_family_new_atleast += 1
	plt.title(newseq_test_name + "\n" + "New Seq has a Signal Greater Than/Equal to " + str(num_family_new_atleast) + " Family Seqs")

	plt.tight_layout()
	plt.savefig(histo_out_dir + newseq_test_name + ".png", bbox_inches = 'tight')

	return


def transform_pair_to_guid(pos_1, pos_1_nt, pos_2, pos_2_nt):
	"""
	"""
	guid = str(pos_1) + "_" + pos_1_nt + "|" + str(pos_2) + pos_2_nt

	return guid


def count_num_edges_per_seq(msa_seqids, msa_seqs, msa_seq_len, \
								observed_edges_w_weights, \
								characters_of_interest, class_max_p_value):
	"""
	"""
	edge_guid_dict = set()
	for pos_1, pos_1_j, pos_2, pos_2_j, pvalue, cont_table in observed_edges_w_weights:
		pos_1_nt = characters_of_interest[pos_1_j]
		pos_2_nt = characters_of_interest[pos_2_j]

		if pvalue > class_max_p_value:
			continue

		guid = transform_pair_to_guid(pos_1, pos_1_nt, pos_2, pos_2_nt)
		edge_guid_dict.add(guid)

	# now, for each sequence, count up all significant edges
	family_seq_tracker = []
	new_seq_tracker = ""
	for seq_index in range(len(msa_seqids)):
		seqid = msa_seqids[seq_index]
		seq = msa_seqs[seq_index]
		num_nucleotides = 0
		num_sig_edges = 0

		for pos_1 in range(len(seq)):
			pos_1_nt = seq[pos_1]

			if pos_1_nt == "-":
				continue

			# else, it's a nucleotide
			num_nucleotides += 1

			# compare to proceeding positions
			for other_pos in range(1, msa_seq_len - pos_1):
				pos_2 = other_pos + pos_1
				pos_2_nt = seq[pos_2]

				if pos_2_nt == "-":
					continue

				guid = transform_pair_to_guid(pos_1, pos_1_nt, pos_2, pos_2_nt)

				if guid in edge_guid_dict:
					num_sig_edges += 1

			# compare to preceeding positions
			for other_pos in range(pos_1):
				pos_2 = other_pos
				pos_2_nt = seq[pos_2]

				if pos_2_nt == "-":
					continue

				guid = transform_pair_to_guid(pos_1, pos_1_nt, pos_2, pos_2_nt)

				if guid in edge_guid_dict:
					num_sig_edges += 1


		# normalize by dividing num_sig_edges by num_nucleotides
		signal = num_sig_edges / num_nucleotides

		if "new" in seqid:
			new_seq_tracker = signal
			print("--------------NEW", signal)
		else:
			family_seq_tracker.append(signal)
			print(seqid, signal)

	return new_seq_tracker, family_seq_tracker
