@testset "Testing matrix_constructors" begin
    @testset "Testing apply_S_plus" begin

    v = ED_sectors.VEC(3.1,5) #5 = [1,0,1]
    w = ED_sectors.VEC(3.1,6) #6 = [0,1,1]
    #println("Testing S plus")
    # println(v)
    # println(apply_S_plus(v,1))
    @test ED_sectors.apply_S_plus(v,1) == ED_sectors.VEC(6.2,7)
    # println(w)
    @test ED_sectors.apply_S_plus(w,1) == ED_sectors.VEC(0, 0)
    @test ED_sectors.is_Null_VEC(ED_sectors.apply_S_plus(w,1))

    end	

    @testset "Testing apply_S_minus" begin
	    v = ED_sectors.VEC(3.1,5) #5 = [1,0,1]
		w = ED_sectors.VEC(3.1,6) #6 = [0,1,1]	    

	    @test ED_sectors.is_Null_VEC(ED_sectors.apply_S_minus(v,1))
	    @test ED_sectors.apply_S_minus(w,1) == ED_sectors.VEC(6.2,4)
    end

    @testset "Testing apply_S_z" begin
	    v = ED_sectors.VEC(3.1,5) #5 = [1,0,1]
		w = ED_sectors.VEC(3.1,6) #6 = [0,1,1]

		@test isapprox(ED_sectors.apply_S_z(v,1).factor, ED_sectors.VEC(-3.1,5).factor)
		@test ED_sectors.apply_S_z(v,1).state==ED_sectors.VEC(-3.1,5).state
		@test w == ED_sectors.VEC(3.1,6)	
        
    end

    @testset "Testing apply_S_x" begin
        v = ED_sectors.VEC(3.1,5) #5 = [1,0,1]
		w = ED_sectors.VEC(3.1,6) #6 = [0,1,1]

	    @test ED_sectors.apply_S_x(v,1) == ED_sectors.VEC(3.1,7)
	    @test ED_sectors.apply_S_x(w,1) == ED_sectors.VEC(3.1,4)
    end

    @testset "Testing apply_S_y" begin
        v = ED_sectors.VEC(3.1,5) #5 = [1,0,1]
		w = ED_sectors.VEC(3.1,6) #6 = [0,1,1]

		@test ED_sectors.apply_S_y(v,1) == ED_sectors.VEC(0.0+3.1*im,7)
		@test ED_sectors.apply_S_y(w,1) == ED_sectors.VEC(-0.0-3.1*im,4)
    end

    @testset "Testing apply_operators" begin
        ops = [ED_sectors.SOP("X",1),ED_sectors.SOP("X",2)]
	    ops4 = [ED_sectors.SOP("-",5)]
	    v = ED_sectors.VEC(1.0,15) #15 = [1,1,1,1]
	    @test ED_sectors.apply_operators(v,ED_sectors.TERM_INTERNAL(2.0,ops)) == ED_sectors.VEC(2.0,9) #9 = [1,0,0,1]
	    @test ED_sectors.is_Null_VEC(ED_sectors.apply_operators(v,ED_sectors.TERM_INTERNAL(3.0,ops4)))
    end

    @testset "Testing agreement of full and symmetry constructors" begin
    	#this will catch a few errors if the full and symmetric matrices differ for some reason

    	#we need some hamiltonian to work with here
        function symbolic_XXZ_Hamiltonian(L :: Int, Delta :: Float64)
		    H = ED_sectors.ABSTRACT_OP(L; name = "XXZ", pbc= true)
		    #H = sum_i sigma_{+}^i sigma_{-}^{i+1}
		    H += ED_sectors.TERM(0.5, "+-")
		    #H = sum_i sigma_{-}^i sigma_{+}^{i+1}
		    H += ED_sectors.TERM(0.5, "-+")
		    #H = sum_i Delta*sigma_z^i sigma_z^{i+1}
		    H += ED_sectors.TERM(Delta*1.0, "ZZ")
		    return H
		end
		L = 2
		abstract_hamiltonian = symbolic_XXZ_Hamiltonian(L,0.8)
		basis = ED_sectors.BASIS(L)
		H1 = ED_sectors.Operator(abstract_hamiltonian,basis)
		println("H1:")
		println(H1)
		H2 = ED_sectors.construct_matrix_full(basis,abstract_hamiltonian)
		println("H2:")
		println(H2)
		@test isapprox(H1,H2)
	end

end