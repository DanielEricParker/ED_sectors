using Nullables

@testset "Testing easy_basis.jl" begin
	@testset "Testing measure_N_up" begin
		# println("Testing N_up")
		f = convert(UInt64, 0) 	#0
		a = convert(UInt64, 1) 	#1
		b = convert(UInt64, 3)	#11
		c = convert(UInt64, 7)	#111
		d = convert(UInt64, 15) #1111
		# display_states([f,a,b,c,d])
		counts = map(s -> ED_sectors.measure_N_up(s,4),[f,a,b,c,d])
		@test counts == [0,1,2,3,4]
	end


	@testset "Testing measure_N_up_A" begin
	    # println("Testing N_up_A")
	    f = convert(UInt64, 0)
	    a = convert(UInt64, 1)
	    b = convert(UInt64, 3)
	    c = convert(UInt64, 7)
	    d = convert(UInt64, 15)
	    # display_states([f,a,b,c,d])
	    counts = map(s -> ED_sectors.measure_N_up_A(s,4),[f,a,b,c,d])
	    @test counts == [0,1,1,2,2]
	end

	@testset "Testing measure_N_up_B" begin
	    # println("Testing N_up_A")
	    f = convert(UInt64, 0)
	    a = convert(UInt64, 1)
	    b = convert(UInt64, 3)
	    c = convert(UInt64, 7)
	    d = convert(UInt64, 15)
	    # display_states([f,a,b,c,d])
	    counts = map(s -> ED_sectors.measure_N_up_B(s,4),[f,a,b,c,d])
	    @test counts == [0,0,1,1,2]
	end

	@testset "Testing make_translation_function." begin
		xx = UInt64(3)
		pf = Complex(1.0)

		tr_fcn = ED_sectors.make_translation_function(10,2,1)
		# println(typeof(tr_fcn))
		#println(digits(x,base=2,pad=10),", ", pf)
		(gx,pf2) = tr_fcn(xx,pf)
		@test digits(gx,base=2,pad=10) == [0,0,1,1,0,0,0,0,0,0]
		@test pf2 ≈ exp(- 2 * pi * im * 1/5)

		tr_fcn = ED_sectors.make_translation_function(10,1,1)
		xx = UInt64(3)
		pf = Complex(1.0)
		(gx,pf2) = tr_fcn(xx,pf)
		@test digits(gx,base=2,pad=10) == [0,1,1,0,0,0,0,0,0,0]
		@test pf2 ≈ exp(- 2 * pi * im * 1/10)
	end	

	@testset "Testing make_Z2B_function" begin
		Z2B_fcn = ED_sectors.make_Z2B_function(5,1)
		Z2B_fcn_2 = ED_sectors.make_Z2B_function(5,-1)
		xxx = UInt64(0)
		pf0 = Complex(1.0)
		# println("x:\t",digits(xxx,base=2,pad=5),", ", pf)

		(gx,pf1) = Z2B_fcn(xxx,pf0)
		@test digits(gx,base=2,pad=5) == [0,1,0,1,0]
		@test pf1 ≈ 1.0 + 0.0*im

		(gx2,pf2) = Z2B_fcn_2(xxx,pf0)
		@test digits(gx2,base=2,pad=5) == [0,1,0,1,0]
		@test pf2 ≈ -1.0 + 0.0*im
	end

	@testset "Testing make_Z2A_function" begin
		Z2A_fcn = ED_sectors.make_Z2A_function(5,1)
		Z2A_fcn_2 = ED_sectors.make_Z2A_function(5,-1)
		xxx = UInt64(0)
		pf0 = Complex(1.0)
		# println("x:\t",digits(xxx,base=2,pad=5),", ", pf)

		(gx,pf1) = Z2A_fcn(xxx,pf0)
		@test digits(gx,base=2,pad=5) == [1,0,1,0,1]
		@test pf1 ≈ 1.0 + 0.0*im

		(gx2,pf2) = Z2A_fcn_2(xxx,pf0)
		@test digits(gx2,base=2,pad=5) == [1,0,1,0,1]
		@test pf2 ≈ -1.0 + 0.0*im
	end

	@testset "Testing make_spin_flip_function" begin
		Z2_fcn = ED_sectors.make_spin_flip_function(5,1)
		Z2_fcn_2 = ED_sectors.make_spin_flip_function(5,-1)
		xxx = UInt64(0)
		xx = UInt64(14)		#[0,1,1,1,0]
		# println(digits(xx,base=2,pad=5))
		pf0 = Complex(1.0)

		(gx,pf1) = Z2_fcn(xxx,pf0)
		@test digits(gx,base=2,pad=5) == [1,1,1,1,1]
		@test pf1 ≈ 1.0 + 0.0*im

		(gx2,pf2) = Z2_fcn_2(xxx,pf0)
		@test digits(gx2,base=2,pad=5) == [1,1,1,1,1]
		@test pf2 ≈ -1.0 + 0.0*im

		(gx3,pf3) = Z2_fcn(xx,pf0)
		@test digits(gx3,base=2,pad=5) == [1,0,0,0,1]
		@test pf3 ≈ 1.0 + 0.0*im

		(gx4,pf4) = Z2_fcn_2(xx,pf0)
		@test digits(gx4,base=2,pad=5) == [1,0,0,0,1]
		@test pf4 ≈ -1.0 + 0.0*im
	end


	@testset "Testing make_Inversion_function" begin
		L = 10
		xx = UInt64(137)
		# println(digits(x,base=2,pad=L)) #[1, 0, 0, 1, 0, 0, 0, 1, 0, 0]
		pf0 = Complex(1.0)

		inv_fcn = ED_sectors.make_Inversion_function(L,1,2)
		inv_fcn_2 = ED_sectors.make_Inversion_function(L,-1,2)

		
		(gx1,pf1) = inv_fcn(xx,pf0)
		@test digits(gx1,base=2,pad=L) == [0,0,0,1,0,0,0,1,1,0]
		@test pf1≈ 1.0 + 0.0*im

		(gx2,pf2) = inv_fcn_2(xx,pf0)
		@test digits(gx2,base=2,pad=L) == [0,0,0,1,0,0,0,1,1,0]
		@test pf2 ≈ -1.0 + 0.0*im

		inv_fcn = ED_sectors.make_Inversion_function(L,1,1)
		inv_fcn_2 = ED_sectors.make_Inversion_function(L,-1,1)

		
		(gx1,pf1) = inv_fcn(xx,pf0)
		@test digits(gx1,base=2,pad=L) == [0,0,01,0,0,0,1,0,0,1]
		@test pf1≈ 1.0 + 0.0*im

		(gx2,pf2) = inv_fcn_2(xx,pf0)
		@test digits(gx2,base=2,pad=L) == [0,0,01,0,0,0,1,0,0,1]
		@test pf2 ≈ -1.0 + 0.0*im
	end



	@testset "Testing make_easy_orbit_function." begin

		# println("\nTesting translation...")
		L = 8 
		# println(symmetries)

		symmetries = Dict("Tr" => 1)
		orb_fcn = ED_sectors.make_easy_orbit_function(L,2,symmetries)
		
		xx = UInt64(1)
		orb = orb_fcn(xx)

		@test orb isa Nullable{ED_sectors.ORBIT}
		if ~isnull(orb)
			elem = get(orb).elements
			for n in 1:length(elem)
				# println(digits(UInt64(2^((2*n)-2)),base=2,pad=L))
				@test elem[UInt64(2^((2*n)-2))] == (-1.0im)^(n-1)
			end
		end


		# println("\nTesting Z2B...")
		symmetries = Dict("PB" => 1)
		orb_fcn = ED_sectors.make_easy_orbit_function(L,2,symmetries)
		x = UInt64(134)
		# println(digits(x,base=2,pad=L)) #[0, 1, 1, 0, 0, 0, 0, 1]
		orb = orb_fcn(x)
		if ~isnull(orb)
			# println("cc_x:", bin(get(orb).representative))
			elem = get(orb).elements
			# println(typeof(elem))
			# for (k,v) in elem
			# 	println(digits(k,base=2,pad=L), ", ", v)
			# end

			keys_digits = sort(map(k -> digits(k,base=2,pad=L),collect(keys(elem))))
			digits_expected = sort([[0, 1, 1, 0, 0, 0, 0, 1],[0, 0, 1, 1, 0, 1, 0, 0]])
			@test keys_digits == digits_expected 

		end

		# println("\nTesting Z2B + translation...")
		
		# symmetries = Dict("Z2B" => -1, "K" => 0)
		# orb_fcn = make_easy_orbit_function(8,2,symmetries)
		

		# x = UInt64(5)
		# orb = orb_fcn(x)

		# println(orb)
		# if ~isnull(orb)
		# 	elem = get(orb).elements
		# 	for (k,v) in elem
		# 		println(bin(k), ", ", v)
		# 	end
		# end



		# println("\nTesting Z2B + Z2A +  translation...")
		L=12
		symmetries = Dict("PB" => -1, "PA" => -1, "Tr" => 1)
		orb_fcn = ED_sectors.make_easy_orbit_function(L,2,symmetries)
		

		x = UInt64(5)
		orb = orb_fcn(x)

		# println(orb)
		if ~isnull(orb)
			elem = get(orb).elements
			# for (k,v) in elem
			# 	println(digits(k,base=2,pad=L), ", ", v)
			# end
			@test length(elem) == 2*2*(12/2)
		end
	end

	#I'm not sure how to stringently test this,
	#but this will at least find at least a few issues I hope

	#I should note that I verified this code pretty carefully against
	#older version, full ED, matlab code other people wrote, etc
	@testset "Testing easy basis function" begin
	    L = 4
	    SzTrZ2basis = ED_sectors.make_basis(L;syms=Dict("SzA" => 2,"Tr" => 0),unitCellSize=2)
	    #println(collect(Set(values(SzTrZ2basis.get_conj_class))))
	    #println("Basis size:", length(SzTrZ2basis.conj_classes))
	    @test length(SzTrZ2basis.conj_classes) == 3
	    # println(SzTrZ2basis.conj_classes)
	    # display_states(collect(keys(SzTrZ2basis.conj_classes)))
	    @test collect(keys(SzTrZ2basis.conj_classes)) == [UInt64(7),UInt64(5),UInt64(15)]
	    # println(SzTrZ2basis.q_numbers)
	    @test SzTrZ2basis.q_numbers == Dict("Tr" => 0,"SzA" => 2)

	    @test_throws AssertionError ED_sectors.make_basis(L;syms=Dict("SzA" => 2,"fdafdsa" => 0),unitCellSize=2)

	end


end