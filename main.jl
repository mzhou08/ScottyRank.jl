using DelimitedFiles
using LinearAlgebra

struct Vertex
  index::UInt32
  name::String
  out_neighbors::Vector{UInt32}
end

n, V = open("graphs/food.txt") do txt
  n = parse(UInt32, readline(txt))
  V = Array{Vertex}(undef, 0)
  for i in 1:n
    data = readdlm(IOBuffer(readline(txt)))
    out_neighbors = Array{UInt32}(undef, 0)
    for j in 2:length(data)
      push!(out_neighbors, convert(UInt32, data[j]))
    end
    push!(V, Vertex(i, data[1], out_neighbors))
  end
  n, V
end

A = zeros(Float64, (n, n))
for i in 1:n
  m = length(V[i].out_neighbors)
  if m == 0
    for j in 1:n
      A[j, i] = 1 / (n - 1);
    end
    A[i, i] = 0
  else
    for j in V[i].out_neighbors
      A[j, i] = 1 / m;
    end
  end
end

display(A)
println();

# damping factor
d = 0.85
M = map(x -> d * x + (1 - d) / n, A)

display(M)
println();
