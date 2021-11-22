using DelimitedFiles
using LinearAlgebra

struct Vertex
  index::UInt32
  name::String
  in_neighbors::Vector{UInt32}
  out_neighbors::Vector{UInt32}
end

# read graph

n, V = open("graphs/food.txt") do txt
  n = parse(UInt32, readline(txt))
  V = Array{Vertex}(undef, 0)
  for i in 1:n
    data = readdlm(IOBuffer(readline(txt)))
    in_neighbors = Array{UInt32}(undef, 0)
    out_neighbors = Array{UInt32}(undef, 0)
    for j in 2:length(data)
      push!(out_neighbors, convert(UInt32, data[j]))
    end
    push!(V, Vertex(i, data[1], in_neighbors, out_neighbors))
  end
  n, V
end

for i in 1:n
  for j in V[i].out_neighbors
    push!(V[j].in_neighbors, i)
  end
end

# PageRank

M = zeros(Float64, (n, n))
for i in 1:n
  m = length(V[i].out_neighbors)
  if m == 0
    for j in 1:n
      M[j, i] = 1 / (n - 1);
    end
    M[i, i] = 0
  else
    for j in V[i].out_neighbors
      M[j, i] = 1 / m;
    end
  end
end

display(M)
println();

# damping factor
d = 0.85
M = map(x -> d * x + (1 - d) / n, M)

display(M)
println();

# HITS

A = zeros(Float64, (n, n))
H = zeros(Float64, (n, n))

for i in 1:n
  for j in V[i].in_neighbors
    A[i, j] = 1
  end
  for j in V[i].out_neighbors
    H[i, j] = 1
  end
end

display(A)
println()

display(H)
println()

a = ones(Float64, n)
h = ones(Float64, n)

function update(a, h) # HITS update function
  na = A * h
  nh = H * a
  na / norm(na), nh / norm(nh)
end

display(a)
println()

display(h)
println()
