module ScottyRank

using DelimitedFiles
using LinearAlgebra
using Printf

export Vertex, Graph
export read_graph
export pagerank_print, pagerank
export hits_print, hits
export generate_adjacency_matrix, generate_adjacency_list

"""
    Vertex

ScottyRank vertex

# Fields
- `index::UInt32`: stores the 1-based index as an unsigned integer
- `in_neighbors::Vector{UInt32}`: stores the indices of incoming neighbors as a list
- `out_neighbors::Vector{UInt32}`: stores the indices of outgoing neighbors as a list
"""
struct Vertex
  index::UInt32
  in_neighbors::Vector{UInt32}
  out_neighbors::Vector{UInt32}
end

"""
    Graph

ScottyRank graph

# Fields
- `num_vertices`: stores the number of vertices as an unsigned integer
- `vertices::Vector{Vertex}`: stores the vertices as a sorted list
"""
struct Graph
  num_vertices::UInt32
  vertices::Vector{Vertex}
end

"""
    read_graph(filepath::String="data/medium-el.txt"); filetype::String="el", zero_index::Bool=false) -> Graph

Reads a graph from an edge list/adjacency list file

# Arguments
- `filepath::String="data/medium-el.txt"`: the path to the source file

# Keywords
- `filetype::String="el"`: "el" for edge list, "al" for adjacency list
- `zero_index::Bool=false`: whether the input file is zero-based

# Returns
- `Graph`: the graph from the source file
"""
function read_graph(filepath::String="data/medium-el.txt";
    filetype::String="el", zero_index::Bool=false)
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
    vertices[i] = Vertex(i, Array{UInt32}(undef, 0), Array{UInt32}(undef, 0))
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
    vertices[i] = Vertex(i, Array{UInt32}(undef, 0), Array{UInt32}(undef, 0))
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

"""
    function pagerank_print(graph::Graph, pg::Vector{Float64}; num_lines::Union{Int64, UInt32}=10, params::Vector{String}=String["vall", "index", "in", "out"]) -> Nothing

Pretty-prints information about the vertices with top PageRank scores to stdout

# Arguments
- `graph::Graph`: the graph
- `pg::Vector{Float64}`: the PageRank scores for the graph

# Keywords
- `num_lines::Union{Int64, UInt32}=10`: the number of vertices whose information is printed
- `params::Vector{String}=String["vall", "index", "in", "out"]`: the types of information printed in order
  - `index`: one-based index
  - `0ndex`: zero-based index
  - `val`: PageRank score, two digits after decimal
  - `vall`: PageRank score, four digits after decimal
  - `valll`: PageRank score, six digits after decimal
  - `in`: number of incoming neighbors
  - `out`: number of outgoing neighbors

# Returns
- `Nothing`
"""
function pagerank_print(graph::Graph, pg::Vector{Float64};
    num_lines::Union{Int64, UInt32}=10, params::Vector{String}=String["vall", "index", "in", "out"])
  if num_lines > graph.num_vertices
    error("invalid num_lines")
  end
  perm = sortperm(pg, rev=true)
  for param in params
    @printf(" - - %5s |", param)
  end
  println()
  println("--- pagerank ---")
  for line in 1:num_lines
    for param in params
      if param == "index"
        @printf("%10d |", perm[line])
      elseif param == "0ndex"
        @printf("%10d |", perm[line] - 1)
      elseif param == "val"
        @printf("%10.2f |", pg[perm[line]])
      elseif param == "vall"
        @printf("%10.4f |", pg[perm[line]])
      elseif param == "valll"
        @printf("%10.6f |", pg[perm[line]])
      elseif param == "in"
        @printf("%10d |", length(graph.vertices[perm[line]].in_neighbors))
      elseif param == "out"
        @printf("%10d |", length(graph.vertices[perm[line]].out_neighbors))
      else
        error("invalid param")
      end
    end
    println()
  end
end

function pagerank(graph::Graph;
    damping::Float64=0.85, modeparam::Tuple{String, Union{Int64, UInt32, Float64}}=("iter", 10))
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
        M[index_to, vertex.index] = 1 / (graph.num_vertices - 1)
      end
      M[vertex.index, vertex.index] = 0
    else
      for index_to in vertex.out_neighbors
        M[index_to, vertex.index] = 1 / num_out_neighbors
      end
    end
  end
  map(x -> damping * x + (1 - damping) / graph.num_vertices, M)
end

"""
    function hits_print(graph::Graph, a::Vector{Float64}, h::Vector{Float64}; num_lines::Union{Int64, UInt32}=10, params::Vector{String}=String["vall", "index", "in", "out"]) -> Nothing

Pretty-prints information about the vertices with top Authority and Hub scores (separately) to stdout

# Arguments
- `graph::Graph`: the graph
- `a::Vector{Float64}`: the Authority scores for the graph
- `h::Vector{Float64}`: the Hub scores for the graph

# Keywords
- `num_lines::Union{Int64, UInt32}=10`: the number of vertices whose information is printed
- `params::Vector{String}=String["vall", "index", "in", "out"]`: the types of information printed in order
  - `index`: one-based index
  - `0ndex`: zero-based index
  - `val`: PageRank score, two digits after decimal
  - `vall`: PageRank score, four digits after decimal
  - `valll`: PageRank score, six digits after decimal
  - `in`: number of incoming neighbors
  - `out`: number of outgoing neighbors

# Returns
- `Nothing`
"""
function hits_print(graph::Graph, a::Vector{Float64}, h::Vector{Float64};
    num_lines::Union{Int64, UInt32}=10, params::Vector{String}=String["vall", "index", "in", "out"])
  if num_lines > graph.num_vertices
    error("invalid num_lines")
  end
  perm_a = sortperm(a, rev=true)
  perm_h = sortperm(h, rev=true)
  for param in params
    @printf(" - - %5s |", param)
  end
  println()
  println("--- authority ---")
  for line in 1:num_lines
    for param in params
      if param == "index"
        @printf("%10d |", perm_a[line])
      elseif param == "0ndex"
        @printf("%10d |", perm_a[line] - 1)
      elseif param == "val"
        @printf("%10.2f |", a[perm_a[line]])
      elseif param == "vall"
        @printf("%10.4f |", a[perm_a[line]])
      elseif param == "valll"
        @printf("%10.6f |", a[perm_a[line]])
      elseif param == "in"
        @printf("%10d |", length(graph.vertices[perm_a[line]].in_neighbors))
      elseif param == "out"
        @printf("%10d |", length(graph.vertices[perm_a[line]].out_neighbors))
      else
        error("invalid param")
      end
    end
    println()
  end
  println("--- hub ---")
  for line in 1:num_lines
    for param in params
      if param == "index"
        @printf("%10d |", perm_h[line])
      elseif param == "0ndex"
        @printf("%10d |", perm_h[line] - 1)
      elseif param == "val"
        @printf("%10.2f |", h[perm_h[line]])
      elseif param == "vall"
        @printf("%10.4f |", h[perm_h[line]])
      elseif param == "valll"
        @printf("%10.6f |", h[perm_h[line]])
      elseif param == "in"
        @printf("%10d |", length(graph.vertices[perm_h[line]].in_neighbors))
      elseif param == "out"
        @printf("%10d |", length(graph.vertices[perm_h[line]].out_neighbors))
      else
        error("invalid param")
      end
    end
    println()
  end
end

function hits(graph::Graph;
    modeparam::Tuple{String, Union{Int64, UInt32, Float64}}=("iter", 10))
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

"""
    function generate_adjacency_matrix(graph::Graph) -> Matrix{Bool}

Generates the adjacency matrix representation for the graph

# Arguments
- `graph::Graph`: the graph

# Returns
- `Matrix{Bool}`: the adjacency matrix representation for the graph
"""
function generate_adjacency_matrix(graph::Graph)
  AM = zeros(Bool, (graph.num_vertices, graph.num_vertices))
  for vertex in graph.vertices, index_to in vertex.out_neighbors
    AM[vertex.index, index_to] = true
  end
  AM
end

"""
    function generate_adjacency_list(graph::Graph) -> Vector{Vector{UInt32}}

Generates the adjacency list representation for the graph

# Arguments
- `graph::Graph`: the graph

# Keywords
- `zero_index::Bool=false`: whether the output file is zero-based

# Returns
- `Vector{Vector{Bool}}`: the adjacency list representation for the graph
"""
function generate_adjacency_list(graph::Graph; zero_index::Bool=false)
  AL = Array{Vector{UInt32}}(undef, graph.num_vertices)
  for vertex in graph.vertices
    AL[vertex.index] = Array{UInt32}(undef, 0)
    for index_to in vertex.out_neighbors
      if zero_index
        push!(AL[vertex.index], index_to - 1)
      else
        push!(AL[vertex.index], index_to)
      end
    end
  end
  AL
end

end
