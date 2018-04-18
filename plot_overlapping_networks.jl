# Given two network files and a threshold, plot the overlapping network

using StatsBase
using MultivariateStats
using NetworkInference
using LightGraphs
using Colors
using GraphPlot
using Compose
using DataFrames
using PyPlot

# ASSUMES each edge only written in one direction
function load_edges(network_file_path, threshold)
	all_edges = readdlm(network_file_path)
	num_edges = Int(round(size(all_edges)[1]*threshold))
	println(num_edges)
	keep_edges = all_edges[1:num_edges, 1:2]
	println(all_edges[1])
	# Make sure all node labels are strings
	for i in 1:length(keep_edges)
		keep_edges[i] = normalize_gene_name(keep_edges[i])
	end
	return num_edges, keep_edges
end

function is_equal_edge(edge1, edge2)
	return (edge1[1] == edge2[1] && edge1[2] == edge2[2]) ||
		(edge1[1] == edge2[2] && edge1[2] == edge2[1])
end

function add_edge!(adjacency_matrix, node_labels_to_indices, edges, i, weight)
	index1 = node_labels_to_indices[edges[i, 1]]
	index2 = node_labels_to_indices[edges[i, 2]]
	adjacency_matrix[index1, index2] += weight
	adjacency_matrix[index2, index1] += weight
end

function normalize_gene_name(gene_name)
    return uppercase(string(gene_name))
end

function get_node_labels(edges1, edges2)
	return sort(unique(vcat(edges1, edges2)))
end

function get_node_labels_to_indices(node_labels, num_nodes)
	return Dict(zip(node_labels, 1:num_nodes))
end

function get_indices_to_node_labels(node_labels, num_nodes)
	return Dict(zip(1:num_nodes, node_labels))
end

function get_group_label_dicts(group_labels_file_path)
	groups = readdlm(group_labels_file_path)
	group_labels = groups[1:1, 1:end]
	num_groups = length(group_labels)
	group_ids = 1:num_groups
    if num_groups == 1
        group_colors = [colorant"lightgreen"]
    else
	    group_colors = linspace(colorant"blue", colorant"red", num_groups)
    end

	group_ids_to_colors = Dict(zip(group_ids, group_colors))

    genes_to_group_ids = Dict{String, Int}()
    for i in 1:num_groups
        for gene in groups[2:end, i:i]
            if gene !== ""
                genes_to_group_ids[normalize_gene_name(gene)] = i
            end
        end 
    end

	return genes_to_group_ids, group_ids_to_colors
end

# Make adjacency matrix
# Add 1 for edges in network1, and 2 for edges in network2
# Stores which edges are in which networks:
# - 0: neither network
# - 1: network 1 only
# - 2: network 2 only
# - 3: both networks
function get_overlap_adjacency_matrix(edges1, num_edges1, edges2, num_edges2, node_labels, num_nodes)
	node_labels_to_indices = get_node_labels_to_indices(node_labels, num_nodes)
	adjacency_matrix = zeros(Int8, (num_nodes, num_nodes))
	for i in 1:num_edges1
		add_edge!(adjacency_matrix, node_labels_to_indices, edges1, i, 1)
	end
	for i in 1:num_edges2
		add_edge!(adjacency_matrix, node_labels_to_indices, edges2, i, 2)
	end
	return adjacency_matrix
end

function get_graph(adjacency_matrix)
	return LightGraphs.SimpleGraphs.SimpleGraph(adjacency_matrix)
end

function get_edge_colors(graph, adjacency_matrix)
	return [
		get_edge_color(adjacency_matrix[edge.src, edge.dst]) for edge in edges(graph)
	]
end

# Sometimes edge is there more than once in each network, in which case weight > 3
# TODO Stop this from happening
# ...But for now, we default weight > 3 to black
function get_edge_color(weight)
	if weight == 1
		return colorant"orange"
	elseif weight == 2
		return colorant"lightblue"
	else
		return colorant"black"
	end
end

function get_communities(graph)
	return label_propagation(graph, 100000)[1]
end

function get_large_communities(communities, min_community_size)
	num_communities = maximum(communities)
	community_sizes = countmap(communities)
	large_communities = []
	for pair in community_sizes
		if pair[2] > min_community_size
			push!(large_communities, pair[1])
		end
	end
	return large_communities
end

function get_community_colors(graph, min_community_size)
	println("Finding communities...")
	communities = get_communities(graph)

	for community in communities
		println(community)
	end

	large_communities = get_large_communities(communities, min_community_size)
	num_large_communities = length(large_communities)

	community_colors = distinguishable_colors(num_large_communities, lchoices = linspace(50, 80, 5), cchoices = linspace(25, 100, 5), hchoices = linspace(0, 340, 10))
	community_colors_dict = Dict(zip(large_communities, community_colors))

	return [community in keys(community_colors_dict) ? community_colors_dict[community] : colorant"white" for community in communities]
end

function get_label_colors(graph, node_labels, num_nodes, group_labels_file_path)
	node_colors = get_default_colors(num_nodes)
	indices_to_node_labels = get_indices_to_node_labels(node_labels, num_nodes)

	genes_to_group_ids, group_ids_to_colors = get_group_label_dicts(group_labels_file_path)
    for i in vertices(graph)
        gene_name = indices_to_node_labels[i]
        if gene_name in keys(genes_to_group_ids)
            group_id = genes_to_group_ids[gene_name]
            node_colors[i] = group_ids_to_colors[group_id]
        end
    end

    return node_colors
end

function get_default_colors(num_nodes)
	return [colorant"white" for i in 1:num_nodes]
end

# color_nodes_by can be:
# - :Community Applies community detection
# - :Group Uses group_labels_file_path
function plot_overlapping_network(network_file_path1, network_file_path2, threshold;
	out_file_path = "", layout = spring_layout, page_size = 150cm,
    node_label_size = 5, node_size = 0.005, group_labels_file_path = "",
    min_community_size = 20, color_nodes_by = :Community,
    locs_x = [], locs_y = [])

	# Get edges for each network
	println("Getting edges for network 1...")
	num_edges1, edges1 = load_edges(network_file_path1, threshold)
	println("Getting edges for network 2...")
	num_edges2, edges2 = load_edges(network_file_path2, threshold)

	plot_overlapping_network(edges1, num_edges1, edges2, num_edges2;
		out_file_path = out_file_path, layout = layout, page_size = page_size,
    	node_label_size = node_label_size, node_size = node_size, group_labels_file_path = group_labels_file_path,
    	min_community_size = min_community_size, color_nodes_by = color_nodes_by,
    	locs_x = locs_x, locs_y = locs_y
    )

    return edges1, num_edges1, edges2, num_edges2
end
function plot_overlapping_network(edges1, num_edges1, edges2, num_edges2;
	out_file_path = "", layout = spring_layout, page_size = 150cm,
    node_label_size = 5, node_size = 0.005, group_labels_file_path = "",
    min_community_size = 20, color_nodes_by = :Community,
    locs_x = [], locs_y = [])

	# Find all nodes labels and give them indices
	println("Finding nodes...")
	node_labels = get_node_labels(edges1, edges2)
	num_nodes = length(node_labels)

	println("Making adjacency matrix...")
	adjacency_matrix = get_overlap_adjacency_matrix(edges1, num_edges1, edges2, num_edges2, node_labels, num_nodes)

	# Get graph and edge colours from adjacency matrix
	println("Making graph...")
	graph = get_graph(adjacency_matrix)
	edge_colors = get_edge_colors(graph, adjacency_matrix)

	# Get node colors
	if color_nodes_by == :Community
		node_colors = get_community_colors(graph, min_community_size)
		for label in node_labels
			println(label)
		end
	elseif color_nodes_by == :Group && group_labels_file_path != ""
		node_colors = get_label_colors(graph, node_labels, num_nodes, group_labels_file_path)
	else
		node_colors = get_default_colors(num_nodes)
	end

	# Get node locations
	if length(locs_x) == 0 && length(locs_y) == 0
		println("Calculating node locations...")
		locs_x, locs_y = layout(graph)
	end

	# Plot the graph
	println("Plotting graph...")
	if length(out_file_path) > 0
		Compose.draw(SVG(out_file_path, page_size, page_size), gplot(graph, locs_x, locs_y, nodelabel = node_labels, nodefillc = node_colors, edgestrokec = edge_colors, NODELABELSIZE = node_label_size, NODESIZE = node_size))
	end

	return locs_x, locs_y
end