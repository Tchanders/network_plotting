# Functions for analysing a network

include("plot_overlapping_networks.jl")

# get_communities(graph) is in plot_overlapping_networks.jl
function get_communities(edges1, num_edges1, edges2, num_edges2)
	graph, node_labels = get_overlap_graph(edges1, num_edges1, edges2, num_edges2)
	communities = get_communities(graph)
	return node_labels, communities
end

function get_nodes_in_community_i(node_labels, communities, i)
	community_node_indices = find(x -> x == i, communities)
	community_node_labels = node_labels[community_node_indices]
	return community_node_labels, community_node_indices
end

function get_nodes_in_community_from_gene(node_labels, communities, gene)
	i = get_community_from_gene(node_labels, communities, gene)
	return get_nodes_in_community_i(node_labels, communities, i)
end

function get_community_from_gene(node_labels, communities, gene)
	return communities[findfirst(node_labels, gene)]
end

function get_overlap_graph(edges1, num_edges1, edges2, num_edges2)
	node_labels = get_node_labels(edges1, edges2)
	num_nodes = length(node_labels)

	adjacency_matrix = get_overlap_adjacency_matrix(edges1, num_edges1, edges2, num_edges2, node_labels, num_nodes)
	graph = get_graph(adjacency_matrix)

	return graph, node_labels
end

function get_community_subgraph(graph, communities, community)
	keep_nodes = find(x -> x == community, communities)
	subgraph, subgraph_ids_to_graph_ids = induced_subgraph(graph, keep_nodes)
	return subgraph, subgraph_ids_to_graph_ids
end

function plot_central_nodes(sub_graph, subgraph_ids_to_graph_ids, node_labels; out_file_path = "", title = "")
	node_degrees = degree(sub_graph)
	node_betweenness = betweenness_centrality(sub_graph)
	fig, ax = plt[:subplots]()

	labels = []
	for i in 1:length(node_degrees)
		push!(labels, node_labels[subgraph_ids_to_graph_ids[i]])
	end

	ax[:scatter](node_degrees, node_betweenness)
	ax[:set_xlabel]("Node degree")
	ax[:set_ylabel]("Betweenness centrality")
	plt[:title](title)
	for (i, label) in enumerate(labels)
		ax[:annotate](label, (node_degrees[i], node_betweenness[i]))
	end
	if length(out_file_path) > 0
		plt[:savefig](out_file_path)
	end
	plt[:close]()
end

function plot_all_central_nodes(edges1, num_edges1, edges2, num_edges2; out_file_pate = "")
	g, nl = get_overlap_graph(edges1, num_edges1, edges2, num_edges2)
	subgraphs = [get_community_subgraph(g, c, i) for i in unique(c)]
	for (i, sg) in enumerate(subgraphs)
    	plot_central_nodes(sg[1], sg[2], nl, out_file_path = "", title = string("Community :", i))
    end
end
