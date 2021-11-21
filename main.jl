using LinearAlgebra;

# Struct for a graph node.
# Index: its unique integer index, ranging from 1 to n.
# Name: the name we assign to the node, e.g. "Wikipedia"
# Neighbors: the nodes that this node points to.
struct node
  index::Int
  name::String
  neighbors::Vector{Int64}
end

open("graphs/food.txt") do io
  n = parse(Int, readline(io))
  for i in 1:n
    line = readline(io)
    data = split(line, " ")
    
    # Initializing the neighbors array
    neighbors = Int[]

    # Constructing the list of neighbors
    # data[i] is the name, so it is not included here
    for j in 2:length(data)
      push!(neighbors, parse(Int, data[j]))
    end

    # Constructing a new node struct
    newNode = node(i,data[1], neighbors)

    println(newNode.index)
    println(newNode.name)
    println(newNode.neighbors)
  end
end
