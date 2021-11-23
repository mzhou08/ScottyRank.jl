using DelimitedFiles
using LinearAlgebra

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

function hits(n::UInt32, V::Vector{Vertex})
  A = zeros(Float64, (n, n))
  H = zeros(Float64, (n, n))
  for i in 1:n
    for j in V[i].out_neighbors
      A[j, i] = 1
    end
    for j in V[i].in_neighbors
      H[j, i] = 1
    end
  end
  A, H
end

function update(a, h) # HITS update function
  na = A * h
  nh = H * a
  na / norm(na), nh / norm(nh)
end

function displayResults(M::Array{Float64, 2}, V::Vector{Vertex}, n::UInt32) 
  for i in 1:5
    M = M*M
  end

  M = M[sortperm(M[:, 1], rev=true), :]

  for i in 1:n
    @printf("%s - %.4f\n", V[i].name, M[i, 1])
  end
end


function main()
  print("Enter file to read, or 'exit' to exit:\n")
  fileName = readline()
  
  while fileName != "exit"
  
    n, V = readgraph(fileName)
    M = pagerank(n, V)
    A, H = hits(n, V)

    # Make sure to update a and h later!
    a = ones(Float64, n)
    h = ones(Float64, n)

    displayResults(M, V, n)

    print("Enter another file to read, or 'exit' to exit:\n")
    fileName = readline()
  end
end

main()
