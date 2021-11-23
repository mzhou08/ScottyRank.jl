module GraphModule

export Vertex, Graph

struct Vertex
  index::UInt32
  num_in_neighbors::UInt32
  in_neighbors::Vector{UInt32}
  num_out_neighbors::UInt32
  out_neighbors::Vector{UInt32}
end

struct Graph
  num_vertices::UInt32
  num_edges::UInt32
  has_names::Bool
  names::Vector{String}
  vertices::Vector{Vertex}
end

end
