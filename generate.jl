using DelimitedFiles
using LinearAlgebra

open("graphs/food.raw") do raw
  open("graphs/food.txt", "w") do txt
    v, e = map(x -> convert(UInt32, x), readdlm(IOBuffer(readline(raw))))

    N = Array{String}(undef, v)
    A = Array{Vector{UInt32}}(undef, v)
    for i in 1:v
      N[i] = readline(raw)
      A[i] = Array{UInt32}(undef, 0)
    end

    for i in 1:e
      t, f = map(x -> convert(UInt32, x), readdlm(IOBuffer(readline(raw))))
      push!(A[t], f)
    end

    write(txt, "$v\n")
    for i in 1:v
      write(txt, "$(N[i])")
      for j in A[i]
        write(txt, " $j")
      end
      write(txt, "\n")
    end
  end
end
