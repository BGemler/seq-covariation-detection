import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np


def generate_color_dict(unique_symbols):
	"""
	this function returns a color map for each unique symbol via index
	"""
	color_map = {}
	colors = cm.rainbow(np.linspace(0, 1, len(unique_symbols)))
	for i in range(len(unique_symbols)):
		color = colors[i]
		color_map[i] = color

	return color_map


def convert_positions_to_circle(num_positions, circle_radius):
	"""
	this function converts a list of positions of length N
	to surround a circle of radius X
	"""
	theta = float(2 * np.pi) / float(num_positions)

	circle_coords = []
	for i in range(num_positions):
		angle = theta * float(i)

		x = float(circle_radius) * math.cos(angle)
		y = float(circle_radius) * math.sin(angle)

		circle_coords.append([x, y])

	return circle_coords


def set_pointcircle_maxradius(cylinder_height, unique_symbols, symbol_points_sizeshrink):
	"""
	this function determines the maximum radius that a symbol point can have
	by limiting the total diameter to fit Y * the number of unqiue symbols
	"""
	cicle_diameter = float(cylinder_height) / float(len(unique_symbols) * symbol_points_sizeshrink)

	return cicle_diameter / 2


def calculate_arrowdelta(initial_x, initial_y, final_x, final_y, initial_z, final_z):
	"""
	Find the dx, dy for the arrow from start, final
	"""
	dx = final_x - initial_x
	dy = final_y - initial_y
	dz = final_z - initial_z

	return dx, dy, dz


def plot_3d_covariation(characters_of_interest, msa_seq_len, \
							max_p_value, very_sig_p_value, \
							circle_radius, cylinder_height, symbol_points_sizeshrink, \
							vector_color, radius_mult, vertical_extension, \
							observed_edges_w_weights, observed_nodes_w_weights):
	"""
	this function plots a series of vectors (e.g., correlation edges) with weights (p-values)
	and next plots initial symbol positions (x, y) with weights (fract)	
	"""
	color_map = generate_color_dict(characters_of_interest)
	circle_coords = convert_positions_to_circle(msa_seq_len, circle_radius)
	max_symbolpoint_radius = set_pointcircle_maxradius(cylinder_height, characters_of_interest, symbol_points_sizeshrink)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection = '3d')

	# add vectors
	num_kept, num_thrown = 0, 0
	num_v_significant = 0
	for position_initial, symbol_index_initial, position_final, symbol_index_final, pvalue, cont_table in observed_edges_w_weights:
		if pvalue > max_p_value:
			num_thrown += 1
			continue

		num_kept += 1

		initial_x, initial_y = circle_coords[position_initial]
		final_x, final_y = circle_coords[position_final]

		initial_z = float(symbol_index_initial) / float(len(characters_of_interest)) * float(cylinder_height)
		final_z = float(symbol_index_final) / float(len(characters_of_interest)) * float(cylinder_height)

		dx, dy, dz = calculate_arrowdelta(initial_x, initial_y, final_x, final_y, initial_z, final_z)

		# if p-value is above very_sig_p_value, set weight to 0.5 (thinner)
		# elif p-value is below very_sig_p_value, set weight to 1.0 (fatter)

		if pvalue > very_sig_p_value:
			vector_weight = 0.25
		else:
			vector_weight = 0.5
			num_v_significant += 1
		
		vector_length = math.sqrt((final_x - initial_x) ** 2 + (final_y - initial_y) ** 2 + (final_z - initial_z) ** 2)
		ax.quiver(initial_x, initial_y, initial_z, dx, dy, dz, linewidths = vector_weight, color = vector_color, \
								arrow_length_ratio = 0)
								#scale_units='xy', scale = 1)


	print("For a p-value threshold of :", max_p_value)
	print("# of vectors kept:", num_kept, " & thrown out", num_thrown)
	print("of vectors kept, ", num_v_significant, " were deemed very significant to p-value threshold:", very_sig_p_value)

	# Add nucleotides
	for position, symbol_index, fract in observed_nodes_w_weights:
		x, y = circle_coords[position]
		z = float(symbol_index) / float(len(characters_of_interest)) * float(cylinder_height)

		radius = float(fract) * max_symbolpoint_radius
		color = color_map[symbol_index]

		ax.scatter(x, y, z, s = radius * radius_mult, c = [color], alpha = 1)

	# Draw each position's vertical, add position text
	for i in range(len(circle_coords)):
		x, y = circle_coords[i]

		initial_x, final_x = x, x
		initial_y, final_y = y, y
		initial_z, final_z = 0, vertical_extension * cylinder_height
		ax.plot(xs = [initial_x, final_x], ys = [initial_y, final_y], zs = [initial_z, final_z], color = "grey")

		position_text = str(i + 1)
		ax.text(initial_x, initial_y, final_z, position_text, color = "yellow")

	# Draw the z-axis labels (e.g., symbols)
	for i in range(len(characters_of_interest)):
		x = min(circle_coords, key = lambda l: l[0])[0]
		y = min(circle_coords, key = lambda l: l[1])[1]

		z = float(i) / float(len(characters_of_interest)) * float(cylinder_height)

		symbol = characters_of_interest[i]
		color = color_map[i]

		ax.text(x, y, z, symbol, color = color, fontsize = 8)


	# graph parameters
	ax.set_facecolor("black")
	ax.grid(False)
	ax.axis('off')
	ax.set_title("Sequence Covariation Anaylsis")

	plt.show()
	plt.draw()

	return
