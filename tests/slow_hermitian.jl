using LinearAlgebra
using SparseArrays
using Arpack

const N = 3
const fill = 10

rand_cols = rand(1:N,fill)
println(rand_cols)
rand_row = [rand(rc:N) for rc in rand_cols]
println(rand_row)

rand_values = rand(Complex{Float64},fill)
println(rand_values)


# for n in 1:length(rand_cols)
# 	println(rand_cols[n],",",rand_row[n])
# end



upper_tri = sparse(rand_cols,rand_row,rand_values)
upper_tri = upper_tri+adjoint(Diagonal(upper_tri))

hermitian_full = upper_tri + adjoint(upper_tri) - Diagonal(upper_tri)
hermitian_view = Hermitian(upper_tri,:U)
@show hermitian_view
println(typeof(hermitian_view))


@show Matrix(hermitian_view) - Matrix(hermitian_full)


