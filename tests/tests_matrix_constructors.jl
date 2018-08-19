include("../main.jl")

using Test


@testset "Testing matrix_constructors" begin
    @testset "Testing apply_S_plus" begin

    v = VEC(3.1,5) #5 = [1,0,1]
    w = VEC(3.1,6) #6 = [0,1,1]
    #println("Testing S plus")
    # println(v)
    # println(apply_S_plus(v,1))
    @test apply_S_plus(v,1) == VEC(6.2,7)
    # println(w)
    @test apply_S_plus(w,1) == VEC(0, 0)
    @test is_Null_VEC(apply_S_plus(w,1))

    end	

    @testset "Testing apply_S_minus" begin
	    v = VEC(3.1,5) #5 = [1,0,1]
		w = VEC(3.1,6) #6 = [0,1,1]	    

	    @test is_Null_VEC(apply_S_minus(v,1))
	    @test apply_S_minus(w,1) == VEC(6.2,4)
    end

    @testset "Testing apply_S_z" begin
	    v = VEC(3.1,5) #5 = [1,0,1]
		w = VEC(3.1,6) #6 = [0,1,1]

		@test isapprox(apply_S_z(v,1).factor, VEC(-3.1,5).factor)
		@test apply_S_z(v,1).state==VEC(-3.1,5).state
		@test w == VEC(3.1,6)	
        
    end

    @testset "Testing apply_S_x" begin
        v = VEC(3.1,5) #5 = [1,0,1]
		w = VEC(3.1,6) #6 = [0,1,1]

	    @test apply_S_x(v,1) == VEC(3.1,7)
	    @test apply_S_x(w,1) == VEC(3.1,4)
    end

    @testset "Testing apply_S_y" begin
        v = VEC(3.1,5) #5 = [1,0,1]
		w = VEC(3.1,6) #6 = [0,1,1]

		@test apply_S_y(v,1) == VEC(0.0+3.1*im,7)
		@test apply_S_y(w,1) == VEC(-0.0-3.1*im,4)
    end

    @testset "Testing apply_operators" begin
        ops = [OP("X",1),OP("X",2)]
	    ops4 = [OP("-",5)]
	    v = VEC(1.0,15) #15 = [1,1,1,1]
	    @test apply_operators(v,TERM(2.0,ops)) == VEC(2.0,9) #9 = [1,0,0,1]
	    @test is_Null_VEC(apply_operators(v,TERM(3.0,ops4)))
    end

    @testset "Testing agreement of full and symmetry constructors" begin
    	#this will catch a few errors if the full and symmetry matrices differ for some reason

    	#we need some hamiltonian to work with here
        function symbolic_XXZ_Hamiltonian(L :: Int, Delta :: Float64)
		    H = ABSTRACT_OP(L, "XXZ", true)
		    #H = sum_i sigma_{+}^i sigma_{-}^{i+1}
		    H += TERMS(0.5, "+-")
		    #H = sum_i sigma_{-}^i sigma_{+}^{i+1}
		    H += TERMS(0.5, "-+")
		    #H = sum_i Delta*sigma_z^i sigma_z^{i+1}
		    H += TERMS(Delta*1.0, "ZZ")
		    return H
		end
		L = 8
		abstract_hamiltonian = symbolic_XXZ_Hamiltonian(L,0.8)
		basis = Basis(L)
		H1 = construct_matrix_sym(basis,abstract_hamiltonian)
		H2 = construct_matrix_full(basis,abstract_hamiltonian)

		@test isapprox(H1,H2)
	end

end