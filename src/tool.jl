# Tool for indexing Hessian
ind_ij(i::Int, j::Int)::Int = (i - 1) * i ÷ 2 + j
ind_nb(nb::Int)::Int = ind_ij(nb + 1, 0)