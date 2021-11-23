module ScottyRank

using DelimitedFiles

export Vertex, Graph, read_edge_list, read_adjacency_list

struct Vertex
  index::UInt32
  in_neighbors::Vector{UInt32}
  out_neighbors::Vector{UInt32}
end

struct Graph
  num_vertices::UInt32
  vertices::Vector{Vertex}
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
  vertices
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
  vertices
end

end
