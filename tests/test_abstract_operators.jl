# H = ED_sectors.ABSTRACT_OP_NEW(L,name="Ising_star",pbc=true)

# #add Z_0 Z_1 + Z_1 Z_2 + Z_2 Z_3 +...
# H += -J*ED_sectors.TERM_NEW["ZZ"]

# #add X_0 + X_1 + ...
# H += -h*ED_sectors.TERM_NEW["X"]

# #add Z_0 X_1 Z_2 + Z_2 X_3 Z_4 + ...
# H += ED_sectors.TERM_NEW[J,"ZXZ",0,repeat=2]

# #add Z_3 X_4 X_5 X_6 Z_7 
# H += ED_sectors.TERM_NEW[3.0,"ZXXXZ",3]



# #S_z^4
# Sz4 = ED_sectors.ABSTRACT_OP_NEW(L,ED_sectors.TERM_NEW["X",4])


####ED_sectors.TERMS_NEW constructor testing
# println(ED_sectors.TERM("ZZ"))

# println(ED_sectors.TERM(3.0,"ZXXXZ",3))

# println(ED_sectors.TERM("ZXZ",5,repeat=2))

# println(ED_sectors.TERM(4,"ZXZ",5,repeat=2))

# println(3.0im*ED_sectors.TERM(4,"ZXZ",5,repeat=2))

# t5 = 3.0im*ED_sectors.TERM(4,"ZXZ",5,repeat=2)
# testop =ED_sectors.ABSTRACT_OP(12, sites="spin-half", name="test op",pbc=false)
# println(ED_sectors.parse_ED_sectors.TERM(testop,t5))

# # t6 = 3.0im*ED_sectors.TERM_NEW(4,"ZXZW",5,repeat=2)
# # ED_sectors.parse_ED_sectors.TERM_new(testop,t6)

# println(testop + t5)


@testset "Testing construction of abstract operators" begin

	@testset "Testing ED_sectors.parse_ED_sectors.TERM" begin
	   @test ED_sectors.parse_term("ZXZIZ",3) == Dict{Int,String}(7=>"Z",4=>"X",3=>"Z",5=>"Z")
	end

    @testset "Testing the TERM constructors" begin
    	t1 = ED_sectors.TERM("ZZ")
    	@test t1.prefactor == ED_sectors.TERM(1, "ZZ", 0; repeat = 1).prefactor
    	@test t1.operators == ED_sectors.TERM(1, "ZZ", 0; repeat = 1).operators
    	@test t1.period == ED_sectors.TERM(1, "ZZ", 0; repeat = 1).period

    	t2 = ED_sectors.TERM(3.0,"ZXXXZ",3)
    	t2b = ED_sectors.TERM(3.0, "ZXXXZ", 3; repeat = 0)

    	@test t2.prefactor == t2b.prefactor
    	@test t2.operators == t2b.operators
    	@test t2.period == t2b.period

    	
		t3 = ED_sectors.TERM("ZXZ",5,repeat=2)
		t3b = ED_sectors.TERM(1, "ZXZ", 5; repeat = 2)
    	@test t3.prefactor == t3b.prefactor
    	@test t3.operators == t3b.operators
    	@test t3.period == t3b.period

		t4 = ED_sectors.TERM(4.3,"ZXZ",5,repeat=2)
		t4b = ED_sectors.TERM(4.3,"ZXZ",5,repeat=2)
		@test t4.prefactor == t4b.prefactor
    	@test t4.operators == t4b.operators
    	@test t4.period == t4b.period
	end

	@testset "Testing scalar multiplication of TERMs" begin
		t4 = 3im*ED_sectors.TERM(4.3,"ZXZ",5,repeat=2)
		t4b = ED_sectors.TERM(12.9im,"ZXZ",5,repeat=2)
		@test isapprox(t4.prefactor, t4b.prefactor)
    	@test t4.operators == t4b.operators
    	@test t4.period == t4b.period

	   
	end

    function is_equal_TERMs(t1 :: ED_sectors.TERM_INTERNAL, t2::ED_sectors.TERM_INTERNAL)
    	if !isapprox(t1.prefactor ,t2.prefactor)
    		println("different prefactor. 1: $(t1.prefactor), 2: $(t2.prefactor)")
    		return false
    	elseif t1.operator != t2.operator
    		println("different operators. 1: $(t1.operator), 2: $(t2.operator)")
    		return false
    	else
    		return true
    	end
    end


    function is_equal_ABSTRACT_OPs(OP1 :: ED_sectors.ABSTRACT_OP, OP2 :: ED_sectors.ABSTRACT_OP)

    	if OP1.name != OP2.name
    		println("unequal names's. 1: $(OP1.name) 2: $(OP2.name)")
    		return false
    	elseif OP1.sites != OP2.sites
    		println("unequal sites's. 1: $(OP1.sites) 2: $(OP2.sites)")
    		return false
    	elseif OP1.pbc != OP2.pbc
    		println("unequal BC's. 1: $(OP1.pbc) 2: $(OP2.pbc)")
    		return false
    	else
    		for n in 1:length(OP1.terms)
    		# if OP1.ED_sectors.TERMs != OP2.ED_sectors.TERMs
    		# println("unequal ED_sectors.TERMs's. 1: $(OP1.ED_sectors.TERMs) 2: $(OP2.ED_sectors.TERMs)")
    			if !is_equal_TERMs(OP1.terms[n],OP2.terms[n])
    				return false
    			end
    		end
    	end
    	return true
    end
	    

	@testset "Testing ABSTRACT_OP" begin
	    # @test ED_sectors.ABSTRACT_OP(10)
	    OP1 = ED_sectors.ABSTRACT_OP(10)
	    OP1b = ED_sectors.ABSTRACT_OP(10, "spin half", "abstract operator", true, Array{ED_sectors.TERM_INTERNAL}([]))
	    @test is_equal_ABSTRACT_OPs(OP1,OP1b)

	    OP2 = ED_sectors.ABSTRACT_OP(10,"X",4)
	    t2 = ED_sectors.TERM_INTERNAL(1.0, [ED_sectors.SOP("X",4)])
	    OP2b = ED_sectors.ABSTRACT_OP(10, "spin half", "abstract operator", true,[t2])
	    @test is_equal_ABSTRACT_OPs(OP2,OP2b)


	    OP3 = ED_sectors.ABSTRACT_OP(10,ED_sectors.TERM("Y",5))
	    t3 = ED_sectors.TERM_INTERNAL(1.0, [ED_sectors.SOP("Y",5)])
		OP3b = ED_sectors.ABSTRACT_OP(10, "spin half", "abstract operator", true,[t3])
		@test is_equal_ABSTRACT_OPs(OP3,OP3b)

		OP4 = ED_sectors.ABSTRACT_OP(10,4.3ED_sectors.TERM("Z",7); name="S_z^7", pbc=false)
		t4 = ED_sectors.TERM_INTERNAL(4.3,[ED_sectors.SOP("Z",7)])
		OP4b = ED_sectors.ABSTRACT_OP(10, "spin half", "abstract operator", true,[t4])

		# H = ED_sectors.ABSTRACT_OP(4; name="Ising Model", pbc=true)
		# H += ED_sectors.TERM("ZZ")
		# H += ED_sectors.TERM("X")
		# println(H)

		# order_parameter = ED_sectors.ABSTRACT_OP(4,0.25*ED_sectors.TERM("ZZ"))
		# println(order_parameter)

	end


	@testset "Testing addition of ABSTRACT_OPs" begin
	    
	    OP1 = ED_sectors.ABSTRACT_OP(10,"X",4)
	    OP2 = ED_sectors.ABSTRACT_OP(10,"Z",5)

	    t3 = [ED_sectors.TERM_INTERNAL(1.0,[ED_sectors.SOP("X",4)]),ED_sectors.TERM_INTERNAL(1.0,[ED_sectors.SOP("Z",5)])]
	    OP3 = ED_sectors.ABSTRACT_OP(10, "spin half", "abstract operator + abstract operator", true,t3)


	    @test is_equal_ABSTRACT_OPs(OP1+OP2, OP3)

	end

	@testset "Testing scalar multiplication of ABSTRACT_OPs" begin
	    OP1 = ED_sectors.ABSTRACT_OP(10,TERM(2,"X",4))
	    lambda = 3
	    lambda_OP1 = ED_sectors.ABSTRACT_OP(10,TERM(6,"X",4))
	    @test is_equal_ABSTRACT_OPs(lambda*OP1,lambda_OP1)
	end


	@testset "Testing addition of ABSTRACT_OPs and TERMs" begin
	    

	    t3 = [ED_sectors.TERM_INTERNAL(1.0,[ED_sectors.SOP("X",4)]),ED_sectors.TERM_INTERNAL(1.0,[ED_sectors.SOP("Z",5)])]
	    OP3 = ED_sectors.ABSTRACT_OP(10, "spin half", "abstract operator + abstract operator", true,t3)


	    t4 = [ED_sectors.TERM_INTERNAL(3.4im,[ED_sectors.SOP("X",4)]),ED_sectors.TERM_INTERNAL(3.4im,[ED_sectors.SOP("Z",5)])]
	    OP4 = ED_sectors.ABSTRACT_OP(10, "spin half", "abstract operator + abstract operator", true,t4)
	    @test is_equal_ABSTRACT_OPs(3.4im*OP3, OP4)

	end


	@testset "Testing ED_sectors.parse_TERM_to_internal" begin
		OP = ED_sectors.ABSTRACT_OP(10)
		# println(OP)

		t2b = ED_sectors.TERM(3, "ZXXZ", 3; repeat = 0)
		@test ED_sectors.parse_TERM_to_internal(OP,t2b)[1].operator == ED_sectors.TERM_INTERNAL(3,[ED_sectors.SOP("Z",3),ED_sectors.SOP("X",4),ED_sectors.SOP("X",5),ED_sectors.SOP("Z",6)]).operator

		t4 = ED_sectors.TERM(4.3,"ZZZ",5,repeat=1)
		@test ED_sectors.parse_TERM_to_internal(OP,t4)[5].operator == [ED_sectors.SOP("Z",0),ED_sectors.SOP("Z",1),ED_sectors.SOP("Z",9)]

		OP2 = ED_sectors.ABSTRACT_OP(10;pbc = false)
		t4 = ED_sectors.TERM(4.3,"ZZZ",5,repeat=1)
		# println(ED_sectors.parse_ED_sectors.TERM_to_internal(OP2,t4))
		@test ED_sectors.parse_TERM_to_internal(OP2,t4)[3].operator == [ED_sectors.SOP("Z",7),ED_sectors.SOP("Z",8),ED_sectors.SOP("Z",9)]


		OP2 = ED_sectors.ABSTRACT_OP(10;pbc = false)
		t4 = ED_sectors.TERM(4.3,Dict(3=>"Z",5=>"Y"),repeat=0)
		# println(ED_sectors.parse_ED_sectors.TERM_to_internal(OP2,t4))
		@test ED_sectors.parse_TERM_to_internal(OP2,t4)[1].operator == [ED_sectors.SOP("Z",3),ED_sectors.SOP("Y",5)]
	end


end