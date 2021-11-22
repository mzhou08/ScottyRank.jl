using DelimitedFiles
using LinearAlgebra

struct Vertex
  index::UInt32
  name::String
  in_neighbors::Vector{UInt32}
  out_neighbors::Vector{UInt32}
end

function readgraph(filename::String)
  file = open("graphs/" * filename * ".txt")
  n = parse(UInt32, readline(file))
  V = Array{Vertex}(undef, 0)
  for i in 1:n
    data = readdlm(IOBuffer(readline(file)))
    in_neighbors = Array{UInt32}(undef, 0)
    out_neighbors = Array{UInt32}(undef, 0)
    for j in 2:length(data)
      push!(out_neighbors, convert(UInt32, data[j]))
    end
    push!(V, Vertex(i, data[1], in_neighbors, out_neighbors))
  end
  for i in 1:n, j in V[i].out_neighbors
    push!(V[j].in_neighbors, i)
  end
  close(file)
  n, V
end

function pagerank(n::UInt32, V::Vector{Vertex})
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
  d = 0.85
  map(x -> d * x + (1 - d) / n, M)
end

n, V = readgraph("food")
M = pagerank(n, V)

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
