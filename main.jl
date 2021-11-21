using LinearAlgebra;

struct node
  index::Int
  name::String
  neighbors::Vector{Int64}
end

open("graphs/test.txt") do io
  n = parse(Int, readline(io))
  for i in 1:n
    line = readline(io)
    data = split(line, " ")
    
    neighbors = Int[]
    for j in 2:length(data)
      push!(neighbors, parse(Int, data[j]))
    end

    newNode = node(i,data[1], neighbors)

    println(newNode.index)
    println(newNode.name)
    println(newNode.neighbors)
  end
end
