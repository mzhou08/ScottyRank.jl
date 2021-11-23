using GraphModule

module TestModule

export test_function

function test_function(V::Vertex)
  V.index, V.num_in_neighbors, V.num_out_neighbors
end

end
