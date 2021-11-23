module ScottyRank

using DelimitedFiles
using LinearAlgebra


export Vertex, Graph

struct Vertex
  index::UInt32
  in_neighbors::Vector{UInt32}
  out_neighbors::Vector{UInt32}
end

struct Graph
  num_vertices::UInt32
  vertices::Vector{Vertex}
end


export read_graph

function read_graph(filepath::String="data/medium-el.txt", filetype::String="el")
  if filetype == "el"
    read_edge_list(filepath)
  elseif filetype == "al"
    read_adjacency_list(filepath)
  else
    error("invalid filetype")
  end
end

function read_edge_list(filepath::String)
  file = open(filepath)
  num_vertices, num_edges = map(x -> convert(UInt32, x), readdlm(IOBuffer(readline(file))))
  vertices = Array{Vertex}(undef, num_vertices)
  for i in 1:num_vertices
    vertices[i] = Vertex(i, Array{UInt32}(undef, 0), Array{UInt32}(undef, 0));
  end
  for _ in 1:num_edges
    index_from, index_to = map(x -> convert(UInt32, x), readdlm(IOBuffer(readline(file))))
    push!(vertices[index_from].out_neighbors, index_to)
    push!(vertices[index_to].in_neighbors, index_from)
  end
  close(file)
  Graph(num_vertices, vertices)
end

function read_adjacency_list(filepath::String)
  file = open(filepath)
  num_vertices = parse(UInt32, readline(file))
  vertices = Array{Vertex}(undef, num_vertices)
  for i in 1:num_vertices
    vertices[i] = Vertex(i, Array{UInt32}(undef, 0), Array{UInt32}(undef, 0));
  end
  for index_from in 1:num_vertices, index_to in map(x -> convert(UInt32, x), readdlm(IOBuffer(readline(file))))
    push!(vertices[index_from].out_neighbors, index_to)
    push!(vertices[index_to].in_neighbors, index_from)
  end
  close(file)
  Graph(num_vertices, vertices)
end


export pagerank

function pagerank(graph::Graph; damping::Float64=0.85, modeparam::Tuple{String, Union{Int64, UInt32, Float64}}=("iter", 10))
  if damping < 0 || damping > 1
    error("invalid damping")
  end
  M = pagerank_matrix(graph, damping)
  if modeparam[1] == "iter"
    if !(isinteger(modeparam[2])) || modeparam[2] < 0
      error("invalid param")
    end
    pagerank_iteration(graph.num_vertices, M, UInt32(modeparam[2]))
  elseif modeparam[1] == "epsi"
    if modeparam[2] <= 0
      error("invalid param")
    end
    pagerank_epsilon(graph.num_vertices, M, Float64(modeparam[2]))
  else
    error("invalid mode")
  end
end

function pagerank_iteration(num_vertices::UInt32, M::Matrix{Float64}, num_iterations::UInt32)
  Base.power_by_squaring(M, num_iterations) * ones(Float64, num_vertices) / num_vertices
end

function pagerank_epsilon(num_vertices::UInt32, M::Matrix{Float64}, epsilon::Float64)
  prev = ones(Float64, num_vertices) / num_vertices
  curr = M * prev
  while norm(prev - curr) > epsilon
    prev, curr = curr, M * curr
  end
  curr
end

function pagerank_matrix(graph::Graph, damping::Float64)
  M = zeros(Float64, (graph.num_vertices, graph.num_vertices))
  for vertex in graph.vertices
    num_out_neighbors = length(vertex.out_neighbors)
    if num_out_neighbors == 0
      for index_to in 1:graph.num_vertices
        M[index_to, vertex.index] = 1 / (graph.num_vertices - 1);
      end
      M[vertex.index, vertex.index] = 0
    else
      for index_to in vertex.out_neighbors
        M[index_to, vertex.index] = 1 / num_out_neighbors;
      end
    end
  end
  map(x -> damping * x + (1 - damping) / graph.num_vertices, M)
end


export hits_matrices, hits_update

function hits_matrices(graph::Graph)
  A = zeros(Float64, (graph.num_vertices, graph.num_vertices))
  H = zeros(Float64, (graph.num_vertices, graph.num_vertices))
  for vertex in graph.vertices
    for index_to in vertex.out_neighbors
      A[index_to, vertex.index] = 1
    end
    for index_from in vertex.in_neighbors
      H[index_from, vertex.index] = 1
    end
  end
  A, H
end

function hits_update(A::Matrix{Float64}, H::Matrix{Float64}, a::Vector{Float64}, h::Vector{Float64})
  normalize(A * h), normalize(H * a)
end

end
