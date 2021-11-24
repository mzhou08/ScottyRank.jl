module ScottyRank

using DelimitedFiles
using LinearAlgebra
using GraphRecipes
using Plots


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

function read_graph(filepath::String="data/medium-el.txt"; filetype::String="el", zero_index::Bool=false)
  if filetype == "el"
    read_edge_list(filepath, zero_index)
  elseif filetype == "al"
    read_adjacency_list(filepath, zero_index)
  else
    error("invalid filetype")
  end
end

function read_edge_list(filepath::String, zero_index::Bool)
  file = open(filepath)
  num_vertices, num_edges = map(x -> convert(UInt32, x), readdlm(IOBuffer(readline(file))))
  vertices = Array{Vertex}(undef, num_vertices)
  for i in 1:num_vertices
    vertices[i] = Vertex(i, Array{UInt32}(undef, 0), Array{UInt32}(undef, 0));
  end
  for _ in 1:num_edges
    index_from, index_to = map(x -> convert(UInt32, x), readdlm(IOBuffer(readline(file))))
    if zero_index
      push!(vertices[index_from + 1].out_neighbors, index_to + 1)
      push!(vertices[index_to + 1].in_neighbors, index_from + 1)
    else
      push!(vertices[index_from].out_neighbors, index_to)
      push!(vertices[index_to].in_neighbors, index_from)
    end
  end
  close(file)
  Graph(num_vertices, vertices)
end

function read_adjacency_list(filepath::String, zero_index::Bool)
  file = open(filepath)
  num_vertices = parse(UInt32, readline(file))
  vertices = Array{Vertex}(undef, num_vertices)
  for i in 1:num_vertices
    vertices[i] = Vertex(i, Array{UInt32}(undef, 0), Array{UInt32}(undef, 0));
  end
  for index_from in 1:num_vertices, index_to in map(x -> convert(UInt32, x), readdlm(IOBuffer(readline(file))))
    if zero_index
      push!(vertices[index_from].out_neighbors, index_to + 1)
      push!(vertices[index_to + 1].in_neighbors, index_from)
    else
      push!(vertices[index_from].out_neighbors, index_to)
      push!(vertices[index_to].in_neighbors, index_from)
    end
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


export hits

function hits(graph::Graph; modeparam::Tuple{String, Union{Int64, UInt32, Float64}}=("iter", 10))
  A, H = hits_matrices(graph)
  if modeparam[1] == "iter"
    if !(isinteger(modeparam[2])) || modeparam[2] < 0
      error("invalid param")
    end
    hits_iteration(graph.num_vertices, A, H, UInt32(modeparam[2]))
  elseif modeparam[1] == "epsi"
    if modeparam[2] <= 0
      error("invalid param")
    end
    hits_epsilon(graph.num_vertices, A, H, Float64(modeparam[2]))
  else
    error("invalid mode")
  end
end

function hits_iteration(num_vertices::UInt32, A::Matrix{Float64}, H::Matrix{Float64}, num_iterations::UInt32)
  a, h = ones(Float64, num_vertices), ones(Float64, num_vertices)
  for _ in 1:num_iterations
    a, h = hits_update(A, H, a, h)
  end
  a, h
end

function hits_epsilon(num_vertices::UInt32, A::Matrix{Float64}, H::Matrix{Float64}, epsilon::Float64)
  prev_a, prev_h = ones(Float64, num_vertices), ones(Float64, num_vertices)
  curr_a, curr_h = hits_update(A, H, prev_a, prev_h)
  while norm(prev_a - curr_a) > epsilon || norm(prev_h - curr_h) > epsilon
    prev_a, prev_h, (curr_a, curr_h) = curr_a, curr_h, hits_update(A, H, curr_a, curr_h)
  end
  curr_a, curr_h
end

function hits_update(A::Matrix{Float64}, H::Matrix{Float64}, a::Vector{Float64}, h::Vector{Float64})
  normalize(A * h), normalize(H * a)
end

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


export visualize

function visualize(graph::Graph)
  AM = zeros(Bool, graph.num_vertices, graph.num_vertices)
  pg = pagerank(graph, modeparam=("iter", 100))
  for vertex in graph.vertices, index_to in vertex.out_neighbors
    AM[vertex.index, index_to] = true
  end
  graphplot(AM,
            method=:chorddiagram,
            names=map(x -> text(string(x), 6), pg),
            linecolor=:black,
            fillcolor=:darkgrey
           )
  #  graphplot(AM,
            #  markersize = 0.2,
            #  node_weights = pg,
            #  markercolor = range(colorant"yellow", stop=colorant"red", length=11),
            #  names = 1:11,
            #  fontsize = 10,
            #  linecolor = :darkgrey
           #  )
end

end
