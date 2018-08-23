var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ED_sectors-Documentation-1",
    "page": "Overview",
    "title": "ED_sectors Documentation",
    "category": "section",
    "text": "Welcome to ED_sectors, a package which provides a fast and straightforward interface for doing exact diagonalization. Its main purpose is to diagonalize Hamiltonians for spin-1/2 chains, taking advantage of symmetries to speed up the computations – hence \"sectors\". ED_sectors is written in Julia, a new language that combines python-like syntax with speed comperable to C and was designed specifically for scientific computation. The code is designed to have easy-to-use syntax while remaining speedy."
},

{
    "location": "index.html#Quick-Example-1",
    "page": "Overview",
    "title": "Quick Example",
    "category": "section",
    "text": "Let\'s start off with a quick example of an ED_sectors program just to see how easy it can be. Let\'s take one of the most popular spin chains, the XXZ model:H = -sum_i=1^L S_i^x S_i+1^x + S_i^y S_i+1^y + Delta S_i^z S_i+1^zThis has several symmetries, including translation and conservation of total spin S_mathrmtot = sum_i S_i^z. We can easily take advantage of those to reduce the dimension of the Hilbert space.L = 16 #length of the spin chain\n\n#create the basis with total Sz = +21and work in the second translation sector\nbasis = BASIS(L; syms = Dict(\"Sz\" => 1, \"Tr\" => 2))\n\n#create an abstract representation of the XXZ Hamiltonian\nXXZ = ABSTRACT_OP(L; name = \"XXZ\", pbc = true)\nXXZ -= TERM(\"XX\") - TERM(\"YY\") - 0.3TERM(\"ZZ\")\n\n#construct an explicit matrix for the operator in the chosen basis\nH = Operator(XXZ, basis)\n\n#get the ground state energy and wavefunction\n(E_0, psi_0) = k_eigvals(H,1)\n\n#compute the magnetization in the ground state\n\nS_tot = Operator(TERM(\"Z\"),basis)\n\nmag_Z = expectation(S_tot,psi_0)In only a few lines we can define a basis, create a Hamiltonian, find its ground state, and measure expectation values. See the Tutorial or the Full API to learn more!"
},

{
    "location": "index.html#Installation-1",
    "page": "Overview",
    "title": "Installation",
    "category": "section",
    "text": "Eventually ED_sectors will be a true Julia package, but for now you need to install it manually. To install ED_Sectors:Install Julia 1.0 or later from https://julialang.org/downloads/. Linux repositories often have older versions of Julia, so it is advised that you use the latest version from the Julia website. Installing Julia is painless: you just drop a folder onto your computer. No linking libraries is necessary!\nClone this package by navigating to a directory of your choice and running  \ngit clone https://github.com/DanielEricParker/ED_sectors\nDownload and install the necessary Julia libraries by running the following commands in the Julia REPL.  \nusing Pkg\nPkg.add(\"Nullables\")\nPkg.add(\"IterativeSolvers\")\nPkg.add(\"Arpack\")\nInclude ED_sectors in your julia program via  \npush!(LOAD_PATH,\"[your path to ED_sectors here]/src/\")\nusing ED_sectors"
},

{
    "location": "index.html#Contents-1",
    "page": "Overview",
    "title": "Contents",
    "category": "section",
    "text": "To learn about the \"moving pieces\" in ED_sectors, please see the Tutorial. More advanced user may wish to consult the Full API documentation. Many common computations are already built-in, including computing entanglement entropy, autocorrelations in time, measurements at temperature, and finding phase transitions! In the case you need to extend the functionality, the internal functions and datatypes are documented under Internals.DocTestSetup = quote\n	using ED_sectors\nendThis package was created by Daniel E. Parker at UC Berkeley. Special thanks to Thomas Scaffidi and Romain Vasseur for inspiring me to create this package and testing out much of its functionality."
},

{
    "location": "tutorial.html#",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial.html#Tutorial-1",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "section",
    "text": "Pages = [\"tutorial.md\"]\nDepth = 2DocTestSetup = quote\n	using ED_sectors\nend"
},

{
    "location": "tutorial.html#Overview,-or-Big-Matrices-are-Hard-1",
    "page": "Tutorial",
    "title": "Overview, or Big Matrices are Hard",
    "category": "section",
    "text": "This tutorial will introduce the basic functionality of the ED_sectors package. It is self-contained from a programming perspective, but assumes some knowledge of quantum mechanics.If you\'re already familiar with spin chains and some of their computational difficulties, you may wish to skip to the Abstract Operators section below.Recall that the time-independent Schrodinger Equation	left H psi_nright = E_n left psi_nrightdetermines the main physical quantities –- the energy eigenvalues and energy eigenstates of our system. Our main task is to solve this eigenvalue problem. For a spin-1/2 chain, the Hilbert space of L is mathcalH = otimes_i=0^L-1 mathbbC^2, whose dimension is N = 2^L. Solving the eigensystem for H naively requires O(N^3) time and O(N^2) space. In practical terms, one cannot solve this beyond about L = 16 on a desktop and perhaps L = 26 on the largest supercomputers, mostly due to the memory constraints. To understand a spin chain, we often want to be as close to the thermodynamic limit, L to infty as possible. Bigger, in this case, is better. However, the infamous exponential increase in the size of the Hilbert space limits us to very small systems. That is, unless we do something more clever."
},

{
    "location": "tutorial.html#How-to-be-clever-1",
    "page": "Tutorial",
    "title": "How to be clever",
    "category": "section",
    "text": "There are several ways we can be more clever and reach larger system sizes. "
},

{
    "location": "tutorial.html#Few-Eigenstate-Approaches-1",
    "page": "Tutorial",
    "title": "Few-Eigenstate Approaches",
    "category": "section",
    "text": "The simplest thing to do is to cut down on what information we need. Oftentimes, the ground state of a system is the only interest part and the rest is irrelevant; computing only a single eigenvalue is much easier than finding 2^L.One can show that, starting with any state leftphiright, then repeatedly acting on it with e^-epsilon H approx 1 - epsilon H for some small epsilon  0, we will eventually converge to the ground state:lim_N to infty left( 1 - epsilon H right)^N leftphiright = left0rightPhysically, multiplication by e^-epsilon H performs imaginary time-evolution, which cools down the system until it is stuck in its ground state. Notice, however, that all we need here is a single vector leftphiright and the matrices we multiply it by are generically very sparse, so they are small to store and quick to multiply. This suggests that we can do much better than the naive O(N^3) expectation above.The Lanczos Algorithm https://en.wikipedia.org/wiki/Lanczos_algorithm does exactly this. It gives a very quick way to find the smallest k eigenvalues and eigenvector while employing only matrix-vector multiplication rather than matrix-matrix multiplication and is therefore much, much faster than finding the full spectrum. (Actual runtime estimates are hard to find, but wikipedia estimates it as O(Nkd) where N is the size of the matrix, k eigenvectors are computed, and d is the density of the matrix.) In practice, on can find the ground states of chains of length L = 22 in a few minutes on a modern desktop.The limit here is usually storing the Hamiltonian matrix rather than the time to perform the computation. Some codes therefore implement \"on-the-fly\" construction of the matrix where the full matrix is never stored, but recomputed one row at a time during each matrix-vector multiplication. This can further extend the number of spins considered."
},

{
    "location": "tutorial.html#Symmetries-1",
    "page": "Tutorial",
    "title": "Symmetries",
    "category": "section",
    "text": "Unfortunately, there are two common situations where this is insufficient. If one cares about dynamical questions, then we often need to multiply by the propagator U(t) = e^-i H t. If we know all the eigenvalues E_n and the (unitary) matrix of eigenvectors V, thenU(t) = V D(t) V^daggerwhere D(t) is the diagonal matrix whose elements are e^-i E_n t. Since multiplying by diagonal matrices is quite fast, this means time-evolution is \"only\" a single matrix-matrix multiplication and, moreover, this is probably about the best one can do. The problem is, finding the full eigensystem is still a very hard operation, so we want to think about how to reduce its complexity as much as possible. Perhaps the best way to do this is to take advantage of symmetries of the problem.note: Note\nAn alternative approach is to give up on exact answers and instead perform approximate time-evolution. It turns out there is a clever way to do this using the Lanczos algorithm as well. Using this, one can reach times around t sim 10^4J where J is the energy scale of a single spin.This is implemented in the quantum dynamite package at https://dynamite.readthedocs.io/en/latest/. This method is particularly amenable to parallelization, enabling systems of size L sim 30.To understand the role of symmetries, let\'s consider an example which will follow us through the rest of the tutorial: the XXZ spin chain:H = -sum_i=1^L S_i^x S_i+1^x + S_i^y S_i+1^y + Delta S_i^z S_i+1^zFor the uninitiated, the XXZ chain is a model for quantum magnetism in one dimension and is one of the most-studied models in condensed matter physics. One of its crucial features is that the total spinS_mathrmtot = sum_i S_i^zis a conserved quantity. (One can check explicitly that S_mathrmtot H = 0.) This symmetry partitions out Hilbert space. The 2^L states each have a representation in the S^z basis with a definite number of spin ups. For instance, the state leftuparrowuparrowdownarrowuparrowdownarrowright has S^z = 3. In the jargon, S^z is a \"good quantum number\" for this model, and has a range -L le S^z le L. There is exactly one state (all up) with S^z = L, L states with S^z = L-1 and, generically, binomLk states with S^z = k, for a total of 2^L.We can then think of our Hilbert space mathcalH as made up of the direct sum of smaller Hilbert spaces mathcalH_k of size binomLk, i.e.mathcalH = mathcalH_-L oplus cdots mathcalH_L-1 oplus mathcalH_LThe fact that S^z is conserved means that the Hamiltonian act on a vector in one S^z-sector mathcalH_a and end up with a vector with any components in some other S^z-sector mathcalH_b. Let\'s think about what this means in terms of explicit matrices. If we choose a basis for mathcalH which is made by concatenating the bases for each of the mathcalH_k in order, then the Hamiltonian will be block diagonal:beginpmatrix\nH_-L 	 0 		 cdots  0\n0		 H_L-1	 cdots  0\nvdots 	 vdots	 ddots  vdots \n0 		 0 	 cdots  0  H_L\nendpmatrixwhere each block matrix H_k is the projection of H to the sector with k spins up.In this form, it is easy to see the benefit: solving the eigenvalue problem for H is now reduced to solving the problem for each of the H_k\'s. For L = 24, the size of the largest block is binom2412 approx 2^21, a substantial reduction in the difficulty of the problem!Here we worked with the particular example of the S_mathrmtot symmetry, but the idea is general. Every (Abelian) symmetry of the problem can be used to further reduce the sizes of the diagonal blocks involved. The most common other symmetries are translation, parity, and inversion. Using some combination of these, one can often reduce the effective size of the system by 2 to 5 spins, and sometimes more. I believe the world record with this method is around L = 48. While a few extra spins might not seem like a lot, consider that most dynamical simulations max out at L = 14. In comparison, working with L = 20 is a huge improvement!The downside of this method is that actually putting the matrix into a block-diagonal form requires a lot of bookkeeping. In fact, the S_mathrmtot example was chosen above since it is particularly simple. Generically, finding the block-diagonal form of the matrix involves a fairly complicated and non-local change of basis. Its form is dictated by the (projective) representations of the finite symmetry groups. note: Note\nThe purpose of this package is to abstract away this troublesome procedure. With ED_sectors, finding a single symmetry block of the Hamiltonian can be done in a single line."
},

{
    "location": "tutorial.html#Abstract-Operators-1",
    "page": "Tutorial",
    "title": "Abstract Operators",
    "category": "section",
    "text": "Let\'s now dispense with the motivation and get down to the business of explaining how this packages works and what its important commands are."
},

{
    "location": "tutorial.html#Choosing-Bases-1",
    "page": "Tutorial",
    "title": "Choosing Bases",
    "category": "section",
    "text": ""
},

{
    "location": "tutorial.html#Making-Matrices-1",
    "page": "Tutorial",
    "title": "Making Matrices",
    "category": "section",
    "text": ""
},

{
    "location": "tutorial.html#Measurements-1",
    "page": "Tutorial",
    "title": "Measurements",
    "category": "section",
    "text": ""
},

{
    "location": "tutorial.html#Dynamics-1",
    "page": "Tutorial",
    "title": "Dynamics",
    "category": "section",
    "text": ""
},

{
    "location": "abstract_operators.html#",
    "page": "Abstract Operators",
    "title": "Abstract Operators",
    "category": "page",
    "text": ""
},

{
    "location": "abstract_operators.html#ED_sectors.TERM",
    "page": "Abstract Operators",
    "title": "ED_sectors.TERM",
    "category": "type",
    "text": "TERM(prefactor,operatorName,firstSite, [repeat=0])\n\nTERM(operatorName) = TERM(1, operatorName, 1, [repeat=1])\nTERM(prefactor,operatorName) = TERM(prefactor, operatorName,1, [repeat=1])\nTERM(operatorName,firstSite,[repeat=0]) = TERM(1, operatorName,1, [repeat=0])\n\nA struct which stores one term in a Hamiltonian. A TERM is an operator whose action maps basis vectors to basis vectors (and not superpositions thereof). For spin-half systems, a term is any string of pauli operators, such as Z_0 Z_1. The TERM struct can be constructed in several different ways fdepending on the situation.\n\nArguments\n\n\'prefactor:: Number\': the numerical prefactor for the term.\n\'operatorName:: String\': a string of the allowed operators. For spin-half, the allowed operators are \"1\", \"I\", \"X\", \"Y\", \"Z\", \"+\", and \"-\".\n\'firstSite:: Int\': the site for the first Pauli operator in the string.\n\'repeat: Int\':: the repeat period for the string of operators. If the firstSite is not specified, then the string repeats with period 1. If firstSite is specified, then the string does not repeat by default –- see examples below.\n\nExamples\n\njulia> TERM(\"ZZ\")\nTERM(1, Dict(0=>\"Z\",1=>\"Z\"), 1)\n\njulia> TERM(3.0,\"ZXXXZ\",3)\nTERM(3.0, Dict(7=>\"Z\",4=>\"X\",3=>\"Z\",5=>\"X\",6=>\"X\"), 0)\n\njulia> TERM(\"ZXZ\",5,repeat=2)\nTERM(1, Dict(7=>\"Z\",5=>\"Z\",6=>\"X\"), 2)\n\njulia> TERM(4.3,\"ZXZ\",5,repeat=2)\nTERM(4.3, Dict(7=>\"Z\",5=>\"Z\",6=>\"X\"), 2)\n\n\n\n\n\n"
},

{
    "location": "abstract_operators.html#Base.:*-Tuple{Number,TERM}",
    "page": "Abstract Operators",
    "title": "Base.:*",
    "category": "method",
    "text": "x * TERM\n*(x,TERM)\n\nScalar multiplication of TERM\'s: the scalar \'x\' acts on the prefactor to \'TERM\'.\n\nExamples\n\njulia> t1 = TERM(4,\"ZXZ\",5,repeat=2)\nTERM(4, Dict(7=>\"Z\",5=>\"Z\",6=>\"X\"), 2)\njulia> 3im*t1\nTERM(0 + 12im, Dict(7=>\"Z\",5=>\"Z\",6=>\"X\"), 2)\n\n\n\n\n\n\n"
},

{
    "location": "abstract_operators.html#ED_sectors.ABSTRACT_OP",
    "page": "Abstract Operators",
    "title": "ED_sectors.ABSTRACT_OP",
    "category": "type",
    "text": "ABSTRACT_OP(L; sites = \"spin half\", name = \"abstract operator\", pbc = true)\nABSTRACT_OP(L, operatorName, site; sites = \"spin half\", name = \"abstract operator\", pbc = true)\nABSTRACT_OP(L, TERM; sites = \"spin half\", name = \"abstract operator\", pbc = true)\n\nThe ABSTRACT_OP struct represents an operator as a string of operator names on sites, along with their numerical prefactors. This also encodes some details of the Hilbert space, including the number of sites L, the type of site in use sites (e.g. spin half), the name for the operator name, and pbc which is true when periodic boundary conditions are employed.\n\nABSTRACT_OP\'s constitute a vector space and can be added and scalar-multiplied by (generically) complex numbers. \n\nExamples\n\nThere are several different types of constructors available for ABSTRACT_OP, enabling the quick definition of \"shells\" of operators which can be used to define complex Hamiltonians, or quick constructors for simple observables.\n\njulia> ABSTRACT_OP(10)\nABSTRACT_OP[name: \"abstract operator\", L: 10, type: spin half, pbc: true, #terms: 0]\n\njulia> ABSTRACT_OP(10,\"X\",4)\nABSTRACT_OP[name: \"abstract operator\", L: 10, type: spin half, pbc: true, #terms: 1]\n1.0 + 0.0im*X_4\n\njulia> ABSTRACT_OP(10,TERM(\"Y\",5))\nABSTRACT_OP[name: \"abstract operator\", L: 10, type: spin half, pbc: true, #terms: 1]\n1.0 + 0.0im*Y_5\n\njulia> ABSTRACT_OP(10,TERM(\"Z\",7); pbc=false)\nABSTRACT_OP[name: \"abstract operator\", L: 10, type: spin half, pbc: false, #terms: 1]\n1.0 + 0.0im*Z_7\n\njulia> ABSTRACT_OP(10,4.3TERM(\"Z\",7); name=\"S_z^7\", pbc=false)\nABSTRACT_OP[name: \"S_z^7\", L: 10, type: spin half, pbc: false, #terms: 1]\n4.3 + 0.0im*Z_7\n\njulia> ABSTRACT_OP(4; name=\"Ising Model\", pbc=true) + TERM(\"ZZ\") + TERM(\"X\")\nABSTRACT_OP[name: \"Ising Model\", L: 4, type: spin half, pbc: true, #terms: 8]\n1.0 + 0.0im*Z_0 Z_1\n1.0 + 0.0im*Z_1 Z_2\n1.0 + 0.0im*Z_2 Z_3\n1.0 + 0.0im*Z_0 Z_3\n1.0 + 0.0im*X_0\n1.0 + 0.0im*X_1\n1.0 + 0.0im*X_2\n1.0 + 0.0im*X_3\n\njulia> order_parameter = ABSTRACT_OP(4,0.25*TERM(\"ZZ\"))\nABSTRACT_OP[name: \"abstract operator\", L: 4, type: spin half, pbc: true, #terms: 4]\n0.25 + 0.0im*Z_0 Z_1\n0.25 + 0.0im*Z_1 Z_2\n0.25 + 0.0im*Z_2 Z_3\n0.25 + 0.0im*Z_0 Z_3\n\nSee also: TERM.\n\n\n\n\n\n"
},

{
    "location": "abstract_operators.html#Base.:*-Tuple{Number,ABSTRACT_OP}",
    "page": "Abstract Operators",
    "title": "Base.:*",
    "category": "method",
    "text": "x * Op\n*(x, Op)\n\nScalar multiplication of ABSTRACT_OPs.\n\nExamples\n\njulia> 	Op = ABSTRACT_OP(10,TERM(2,\"X\",4))\nABSTRACT_OP[name: \"abstract operator\", L: 10, type: spin half, pbc: true, #terms: 1]\n2.0 + 0.0im*X_4\n\njulia> 3*Op\nABSTRACT_OP[name: \"abstract operator\", L: 10, type: spin half, pbc: true, #terms: 1]\n6.0 + 0.0im*X_4\n\n\n\n\n\n\n"
},

{
    "location": "abstract_operators.html#Base.:+-Tuple{ABSTRACT_OP,ABSTRACT_OP}",
    "page": "Abstract Operators",
    "title": "Base.:+",
    "category": "method",
    "text": "O1 + O2\n+(O1,O2)\n\nAdds two ABSTRACT_OPs. Operators must have the same length, site type, and boundary conditions.\n\nExamples\n\njulia> 	OP1 = ABSTRACT_OP(10,\"X\",4); OP2 = ABSTRACT_OP(10,\"Z\",5);\n\njulia> OP1+OP2\nABSTRACT_OP[name: \"abstract operator + abstract operator\", L: 10, type: spin half, pbc: true, #terms: 2]\n1.0 + 0.0im*X_4\n1.0 + 0.0im*Z_5\n\n\n\n\n\n"
},

{
    "location": "abstract_operators.html#Base.:+-Tuple{ABSTRACT_OP,TERM}",
    "page": "Abstract Operators",
    "title": "Base.:+",
    "category": "method",
    "text": "Op + term\n+(Op, term)\n\nAdds a new TERM to an ABSTRACT_OP. The TERM must fit inside the number of sites for the operator.\n\nExamples\n\njulia> ABSTRACT_OP(4; name=\"Ising Model\", pbc=true) + TERM(\"ZZ\") + TERM(\"X\")\nABSTRACT_OP[name: \"Ising Model\", L: 4, type: spin half, pbc: true, #terms: 8]\n1.0 + 0.0im*Z_0 Z_1\n1.0 + 0.0im*Z_1 Z_2\n1.0 + 0.0im*Z_2 Z_3\n1.0 + 0.0im*Z_0 Z_3\n1.0 + 0.0im*X_0\n1.0 + 0.0im*X_1\n1.0 + 0.0im*X_2\n1.0 + 0.0im*X_3\n\n\n\n\n\n"
},

{
    "location": "abstract_operators.html#Abstract-Operators-1",
    "page": "Abstract Operators",
    "title": "Abstract Operators",
    "category": "section",
    "text": "DocTestSetup = quote\n	using ED_sectors\nendTERM\n*(::Number, ::TERM)\nABSTRACT_OP\n*(::Number, ::ABSTRACT_OP)\n+(::ABSTRACT_OP,::ABSTRACT_OP)\n+(:: ABSTRACT_OP, :: TERM)"
},

{
    "location": "basis.html#",
    "page": "Making Bases",
    "title": "Making Bases",
    "category": "page",
    "text": ""
},

{
    "location": "basis.html#ED_sectors.BASIS",
    "page": "Making Bases",
    "title": "ED_sectors.BASIS",
    "category": "type",
    "text": "	BASIS(L;\n	a :: Int64 = 1,\n	syms :: Dict{String,Int} = (),\n	constraint :: Function = identity\n	)\n\nA struct that stores a basis for a Hilbert space. A BASIS is constructed by specifiying its attributes. 	- The number of spins L 	- The number of spins per unit cell a, whose default is 1. 	- Any symmetries of the basis, specified as a dictionary syms. The default is no symmetries. 	- A constraint function state :: UInt64 -> Bool where only states where the constraint evaluates to true are admitted to the Hilbert space.\n\nSymmetries are specified by giving their name and sector. For instance, \"tr\" => 3 represents the translation symmetry in sector 3, where one application of the symmetry moves the state around by 3 basis vectors. The complete list of implemented symmetries is below.\n\nSymmetry Full Name Values\nTr Translation 1 le Tr le La\nP Parity (Spin flip) -11\nI Inversion -11\nSz Total S^z -L le Sz le L\nPA P for even spins -11\nPB P for odd spins -11\nSzA sum_imathrmeven S^z_i -L2 le SzA le L2\nSzB sum_imathrmodd S^z_i -L2 le SzB le L2\n\nUsing some symmetries requires extra conditions. Explicitly, translation symmetry requires periodic boundary conditions, and even/odd parity and spin sectors require a even number of total sites. However, multiple symmetries can be used at the same time.\n\nwarning: Index Convention\nSpins are zero-indexed, so the left-most site is spin zero and is in the \"A\" sublattice.\n\nUsers can also supply arbitrary constraint functions UInt64 -> Bool, which specify states that are allowed in the Hilbert space. \n\nnote: Note\nA BASIS can require a very large amount of memory to store. Internally, it is composed of one hashmap from states in the full Hilbert space to a representative of their conjugacy classes, and another from representatives of conjugacy classes to indices in the small Hilbert space. Altogether, this requires 2^L space. For truly large sizes, one can –- in principle –- generate these hashmaps as needed, but this is not implemented yet here.\n\nExamples\n\njulia> BASIS(8)\nBASIS[L: 8, a: 1, dim: 256, symmetries: none]\n\njulia> BASIS(8; a=2)\nBASIS[L: 8, a: 2, dim: 256, symmetries: none]\n\njulia> BASIS(8; syms = Dict(\"Tr\"=>1))\nBASIS[L: 8, a: 1, dim: 30, symmetries: \"Tr\" => 1]\n\njulia> BASIS(8; a=2, syms = Dict(\"Tr\"=>1))\nBASIS[L: 8, a: 2, dim: 60, symmetries: \"Tr\" => 1]\n\njulia> BASIS(8; syms = Dict(\"Sz\" => 3, \"Tr\" => 2))\nBASIS[L: 8, a: 1, dim: 7, symmetries: \"Sz\" => 3, \"Tr\" => 2]\n\njulia> cons_fcn = x -> x < 17;\n\njulia> BASIS(8; constraint = cons_fcn)\nBASIS[L: 8, a: 1, dim: 17, symmetries: \"constrained\" => 1]\n\n\n\n\n\n"
},

{
    "location": "basis.html#Basis-1",
    "page": "Making Bases",
    "title": "Basis",
    "category": "section",
    "text": "DocTestSetup = quote\n	using ED_sectors\nendBASIS"
},

{
    "location": "matrix_constructors.html#",
    "page": "Operators",
    "title": "Operators",
    "category": "page",
    "text": ""
},

{
    "location": "matrix_constructors.html#ED_sectors.Operator",
    "page": "Operators",
    "title": "ED_sectors.Operator",
    "category": "function",
    "text": "Operator(abstract_op, basis)\nOperator(term, basis)\nOperator(operatorName, site, basis)\n\nConstructs an operator in a given BASIS. One can specify a full ABSTRACT_OP or, as a shortcut for simple obserables, just a single TERM or even the name and site of the operator.\n\nReturns a sparse matrix.\n\nExamples\n\nSimple observables are simple to create.\n\njulia> L = 4; basis = BASIS(L);\n\njulia> Sx2 = Operator(\"X\", 2, basis)\n16×16 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 16 stored entries:\n  [5 ,  1]  =  1.0+0.0im\n  [6 ,  2]  =  1.0+0.0im\n  [7 ,  3]  =  1.0+0.0im\n  [8 ,  4]  =  1.0+0.0im\n  [1 ,  5]  =  1.0+0.0im\n  [2 ,  6]  =  1.0+0.0im\n  [3 ,  7]  =  1.0+0.0im\n  [4 ,  8]  =  1.0+0.0im\n  [13,  9]  =  1.0+0.0im\n  [14, 10]  =  1.0+0.0im\n  [15, 11]  =  1.0+0.0im\n  [16, 12]  =  1.0+0.0im\n  [9 , 13]  =  1.0+0.0im\n  [10, 14]  =  1.0+0.0im\n  [11, 15]  =  1.0+0.0im\n  [12, 16]  =  1.0+0.0im\n\njulia> magnetic_order_parameter = Operator((1/L)*TERM(\"Z\"), basis)\n16×16 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 10 stored entries:\n  [1 ,  1]  =  1.0+0.0im\n  [2 ,  2]  =  0.5+0.0im\n  [3 ,  3]  =  0.5+0.0im\n  [5 ,  5]  =  0.5+0.0im\n  [8 ,  8]  =  -0.5+0.0im\n  [9 ,  9]  =  0.5+0.0im\n  [12, 12]  =  -0.5+0.0im\n  [14, 14]  =  -0.5+0.0im\n  [15, 15]  =  -0.5+0.0im\n  [16, 16]  =  -1.0+0.0im\n\n\nCreating more involved operators, like most Hamiltonians, is also straightforward.\n\njulia> L = 4; basis = BASIS(L);\n\njulia> ising = ABSTRACT_OP(L; name=\"Ising Model\", pbc=true) + 2TERM(\"ZZ\") + TERM(\"X\");\n\njulia> H = Operator(ising,basis)\n16×16 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 68 stored entries:\n  [1 ,  1]  =  8.0+0.0im\n  [2 ,  1]  =  1.0+0.0im\n  [3 ,  1]  =  1.0+0.0im\n  [5 ,  1]  =  1.0+0.0im\n  [9 ,  1]  =  1.0+0.0im\n  [1 ,  2]  =  1.0+0.0im\n  [4 ,  2]  =  1.0+0.0im\n  [6 ,  2]  =  1.0+0.0im\n  [10,  2]  =  1.0+0.0im\n  ⋮\n  [7 , 15]  =  1.0+0.0im\n  [11, 15]  =  1.0+0.0im\n  [13, 15]  =  1.0+0.0im\n  [16, 15]  =  1.0+0.0im\n  [8 , 16]  =  1.0+0.0im\n  [12, 16]  =  1.0+0.0im\n  [14, 16]  =  1.0+0.0im\n  [15, 16]  =  1.0+0.0im\n  [16, 16]  =  8.0+0.0im\n\n\nSee also: ABSTRACT_OP, BASIS.\n\n\n\n\n\n"
},

{
    "location": "matrix_constructors.html#Matrix-Constructors-1",
    "page": "Operators",
    "title": "Matrix Constructors",
    "category": "section",
    "text": "DocTestSetup = quote\n	using ED_sectors\nendOperator"
},

{
    "location": "measurement.html#",
    "page": "Ground State Measurements",
    "title": "Ground State Measurements",
    "category": "page",
    "text": ""
},

{
    "location": "measurement.html#ED_sectors.k_eigvals",
    "page": "Ground State Measurements",
    "title": "ED_sectors.k_eigvals",
    "category": "function",
    "text": "k_eigvals(H,k)\n\nComputes the k lowest-energy eigenvalues for the given Hamiltonian H. \n\nnote: Note\nThis is a simple wrapper for Julia\'s eigs function, which uses the Lanczos algorithm. While this works with either sparse or full matrices for H, sparse matrices are much faster. The cost of the algorithm increases with k and should be only used when k ll operatornamedim(H). If you only need the ground state, then just use k=1.\n\n\n\n\n\n"
},

{
    "location": "measurement.html#ED_sectors.k_eigsys",
    "page": "Ground State Measurements",
    "title": "ED_sectors.k_eigsys",
    "category": "function",
    "text": "k_eigsys(H, k)\n\nComputes the k lowest-energy eigenvalues and eigenvectors for the given Hamiltonian H. \n\n\n\n\n\n"
},

{
    "location": "measurement.html#ED_sectors.expectation",
    "page": "Ground State Measurements",
    "title": "ED_sectors.expectation",
    "category": "function",
    "text": "expectation(psi, O)\n	expectation_rho(rho, O)\n\nFor a state psi, computes the expectation value leftpsiOpsiright  and operator O. This is a simple wrapper for dot(psi, O*psi) and works for any type of matrix O.\n\nFor a density matrix rho, computes the expectation value operatornameTrrho O.\n\n\n\n\n\n"
},

{
    "location": "measurement.html#ED_sectors.reduce_density_matrix",
    "page": "Ground State Measurements",
    "title": "ED_sectors.reduce_density_matrix",
    "category": "function",
    "text": "reduce_density_matrix(phi, l)\n\nComputes the reduced density matrix for the state psi by tracing out all but the first l spins.\n\nwarning: Basis restriction\nAt present this only works when the basis has no symmetries. Support for reduced density matrices with symmetries is coming soon.\n\n\n\n\n\n"
},

{
    "location": "measurement.html#ED_sectors.entanglement_entropy",
    "page": "Ground State Measurements",
    "title": "ED_sectors.entanglement_entropy",
    "category": "function",
    "text": "entanglement_entropy(psi, l; epsilon = 10e-15)\n\nComputes the entanglement entropy S(lL) for a state psi. Removes eigenvalues smaller than a numerical cutoff because they often give underflow errors. \n\nThis uses reduce_density_matrix, so it only works without symmetries.\n\n\n\n\n\n"
},

{
    "location": "measurement.html#ED_sectors.measure_EE",
    "page": "Ground State Measurements",
    "title": "ED_sectors.measure_EE",
    "category": "function",
    "text": "measure_EE(psi, basis)\n\nComputes the entanglement entropy at all cuts in a wavefunction psi in a basis basis. \n\n\n\n\n\n"
},

{
    "location": "measurement.html#Measurement-1",
    "page": "Ground State Measurements",
    "title": "Measurement",
    "category": "section",
    "text": "DocTestSetup = quote\n	using ED_sectors\nendk_eigvals\nk_eigsys\nexpectation\nreduce_density_matrix\nentanglement_entropy\nmeasure_EE"
},

{
    "location": "full_spec.html#",
    "page": "Full Spectrum Measurements",
    "title": "Full Spectrum Measurements",
    "category": "page",
    "text": ""
},

{
    "location": "full_spec.html#ED_sectors.EIGSYS",
    "page": "Full Spectrum Measurements",
    "title": "ED_sectors.EIGSYS",
    "category": "type",
    "text": "EIGSYS(evs, U)\nEIGSYS(H)\n\nA struct for storing the full spectrum evs and eigenstates of a system U. Explicitly, U is the (unitary) matrix which gives the change of basis from the original to energy bases, i.e. the columns of U are the eigenvectors of the system. So for an EIGSYS ES, the full spectrum is ES.evs and the 3rd eigenvector, for example, is ES.U[:,3].\n\nExamples\n\njulia> L = 4; basis = BASIS(L);\n\njulia> ising = ABSTRACT_OP(L; name=\"Ising Model\", pbc=true) + 2TERM(\"ZZ\") + TERM(\"X\");\n\njulia> H = Operator(ising,basis);\n\njulia> ES = EIGSYS(H)\nEIGSYS[dim: 16]\n\njulia> psi3 = ES.U[:,3]\n16-element Array{Complex{Float64},1}:\n  0.08727992975105786 + 0.0im\n -0.23235256582864727 + 0.0im\n  -0.2323525658286471 + 0.0im\n   0.3509044118546869 + 0.0im\n -0.23235256582864688 + 0.0im\n -0.17367654405779392 + 0.0im\n  0.35090441185468646 + 0.0im\n -0.23235256582864808 + 0.0im\n -0.23235256582864716 + 0.0im\n   0.3509044118546871 + 0.0im\n  -0.1736765440577938 + 0.0im\n -0.23235256582864866 + 0.0im\n  0.35090441185468646 + 0.0im\n -0.23235256582864844 + 0.0im\n -0.23235256582864783 + 0.0im\n  0.08727992975105792 - 0.0im\n\njulia> isapprox((H * psi3),(ES.evs[3]*psi3)) #they are the same to numerical error\ntrue\n\n\n\n\n\n\n"
},

{
    "location": "full_spec.html#ED_sectors.thermal_density_matrix",
    "page": "Full Spectrum Measurements",
    "title": "ED_sectors.thermal_density_matrix",
    "category": "function",
    "text": "thermal_density_matrix(eigsys, beta)\n\nComputes the thermal density matrix in the energy basis for a given eigsys and inverse temperature beta. When beta = 0, this gives the \"infinite-temperature\" density matrix.\n\nwarning: Warning\nIt is easy to get underflow or overflow errors by exponentiating large numbers. The largest Float64 is roughly 10^308 sim e^709. It is recommended to double-check that your results make sense after using this function.\n\n\n\n\n\n"
},

{
    "location": "full_spec.html#ED_sectors.r_statistic",
    "page": "Full Spectrum Measurements",
    "title": "ED_sectors.r_statistic",
    "category": "function",
    "text": "r_statistic(eigenvalues :: Array{Float64})\nr_statistic(ES :: EIGSYS)\n\nComputes the R-statistic for a spectrum. The R-statistic is a measure of how close to integrability the system is. For integrable systems, leftrright sim 0386. For GOE level-statistics, leftrright sim 05295. \n\nwarning: Warning\nUsing the R statistic correctly is subtle. It only makes sense systems sufficiently large as to be close to the thermodynamic limit. Moreover, one must restrict to a single symmetry sector. Any residual symmetry will give an artificially low R-statistic. See, e.g., https://arxiv.org/abs/cond-mat/0610854 for details..\n\n\n\n\n\n"
},

{
    "location": "full_spec.html#Full-Spectrum-Measurements-1",
    "page": "Full Spectrum Measurements",
    "title": "Full Spectrum Measurements",
    "category": "section",
    "text": "DocTestSetup = quote\n	using ED_sectors\nendEIGSYSthermal_density_matrix\nr_statistic"
},

{
    "location": "dynamics.html#",
    "page": "Dynamics",
    "title": "Dynamics",
    "category": "page",
    "text": ""
},

{
    "location": "dynamics.html#ED_sectors.evolve_state",
    "page": "Dynamics",
    "title": "ED_sectors.evolve_state",
    "category": "function",
    "text": "evolve_state(psi, eigsys, t)\n\nA function which computes the time-evolution of a state psi under a (already diagonalized) Hamiltonian eigsys for time t. Explicitly, this computes lefte^- i H t psiright.\n\n\n\n\n\n"
},

{
    "location": "dynamics.html#ED_sectors.expectation_time",
    "page": "Dynamics",
    "title": "ED_sectors.expectation_time",
    "category": "function",
    "text": "expectation_time(O, t, psi, eigsys)\nexpectation_time(O, t, rho, eigsys\n\nFor a state psi, measures the expectation of the operator O at time t  evolved under eigsys: leftpsiO(t)psiright. \n\nFor a density operator rho, measures the expectation of the operator O at time t  evolved under eigsys: operatornameTrrho O(t).\n\n\n\n\n\n"
},

{
    "location": "dynamics.html#ED_sectors.autocorrelation",
    "page": "Dynamics",
    "title": "ED_sectors.autocorrelation",
    "category": "function",
    "text": "autocorrelation(O, t, rho, eigsys)\n\nComputes the autocorrelation of an operator O at time t against the density matrix rho when evolved under eigsys: operatornameTrrho O(t) O(0).\n\n\n\n\n\n"
},

{
    "location": "dynamics.html#Dynamics-1",
    "page": "Dynamics",
    "title": "Dynamics",
    "category": "section",
    "text": "DocTestSetup = quote\n	using ED_sectors\nendevolve_state\nexpectation_time\nautocorrelation"
},

{
    "location": "utilities.html#",
    "page": "Utilities",
    "title": "Utilities",
    "category": "page",
    "text": ""
},

{
    "location": "utilities.html#ED_sectors.export_data",
    "page": "Utilities",
    "title": "ED_sectors.export_data",
    "category": "function",
    "text": "export_data(data,name,parameters)\n\nExports arrays to csv files, with parameters in the filename. Particularly useful when varying 1 or more parameters. #Arguments\n\n\'data :: Array{Float64, N} where N\': data to export, an Array\n\'name :: String\': name/prefix for the data, e.g. \"XXZ_eigenvalues\"\n\'parameters :: Array{Tuple{String,Number}}\': Array of Tuples (parameter, value), e.g. (\"Delta\", 0.6)\n\n\n\n\n\n"
},

{
    "location": "utilities.html#ED_sectors.make_filename",
    "page": "Utilities",
    "title": "ED_sectors.make_filename",
    "category": "function",
    "text": "make_filename(name,parameters)\n\nGiven the prefix \'name\' and the parameter/value pairs, generates a filename for the data.\n\n#Arguments\n\n\'name :: String\': name/prefix for the data, e.g. \"XXZ_eigenvalues\"\n\'parameters :: Array{Tuple{String,Number}}\': Array of Tuples (parameter, value), e.g. (\"Delta\", 0.6)\n\n\n\n\n\n"
},

{
    "location": "utilities.html#Utilities-1",
    "page": "Utilities",
    "title": "Utilities",
    "category": "section",
    "text": "DocTestSetup = quote\n	using ED_sectors\nendED_sectors.export_data\nED_sectors.make_filename"
},

{
    "location": "internals.html#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "internals.html#Internals-1",
    "page": "Internals",
    "title": "Internals",
    "category": "section",
    "text": "DocTestSetup = quote\n	using ED_sectors\nend\n"
},

{
    "location": "internals.html#ED_sectors.SOP",
    "page": "Internals",
    "title": "ED_sectors.SOP",
    "category": "type",
    "text": "SOP(name,site)\n\nA struct for storing a single-site operator. This is an internal method which users shouldn\'t have to access.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.TERM_INTERNAL",
    "page": "Internals",
    "title": "ED_sectors.TERM_INTERNAL",
    "category": "type",
    "text": "TERM_INTERNAL(prefactor, operator)\n\nThe internal representation of a term as a (generically complex) prefactor and an array of single-single operators called operator.\n\n\n\n\n\n"
},

{
    "location": "internals.html#Base.:*-Tuple{Number,ED_sectors.TERM_INTERNAL}",
    "page": "Internals",
    "title": "Base.:*",
    "category": "method",
    "text": "*(x,TERM_INTERNAL)\nx*TERM\n\nDefined scalar multiplication of TERM\'s: the scalar acts on the prefactor.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.parse_TERM_to_internal-Tuple{ABSTRACT_OP,TERM}",
    "page": "Internals",
    "title": "ED_sectors.parse_TERM_to_internal",
    "category": "method",
    "text": "parse_term(op, term)\n\nInternal function to parse a TERM into an array of TERM_INTERNALs so they can be added onto an operator.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.parse_term-Tuple{String,Int64}",
    "page": "Internals",
    "title": "ED_sectors.parse_term",
    "category": "method",
    "text": "parse_term(op, term)\n\nInternal function to parse the operators given in the TERM constructor into a Dictionary.\n\n\n\n\n\n"
},

{
    "location": "internals.html#abstract_operators.jl-1",
    "page": "Internals",
    "title": "abstract_operators.jl",
    "category": "section",
    "text": "Modules = [ED_sectors]\nPages = [\"abstract_operators.jl\"]\nPublic = false\nOrder   = [:constant,:type,:function]"
},

{
    "location": "internals.html#ED_sectors.BASISVECTOR",
    "page": "Internals",
    "title": "ED_sectors.BASISVECTOR",
    "category": "type",
    "text": "BASISVECTOR(conj_class, phase_factor)\n\nA struct that links a basis vector to its conjugacy class `conjclass\' and gives the phase factor \'phasefactor\' between them due to the representation.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.CONJCLASS",
    "page": "Internals",
    "title": "ED_sectors.CONJCLASS",
    "category": "type",
    "text": "CONJCLASS(index,norm)\n\nA struct that describes the conjugacy class of a basis element. With symmetries, each conjugacy class is associated to a basis vector. The index is the index of the associated basis vector and the norm is the norm squared of that conjugacy class, so it can be properly normalized.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.ORBIT",
    "page": "Internals",
    "title": "ED_sectors.ORBIT",
    "category": "type",
    "text": "ORBIT(norm, representative, elements)\n\nA struct which stores the orbit of a single element under an Abelian group action.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.check_same_basis_debug-Tuple{Any,Any}",
    "page": "Internals",
    "title": "ED_sectors.check_same_basis_debug",
    "category": "method",
    "text": "check_same_basis_debug(basis1,basis2)\n\nChecks if two bases are the same. UNTESTED! Use for debugging only.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.get_valid_state_function-Tuple{Int64,Int64,Dict{String,Int64}}",
    "page": "Internals",
    "title": "ED_sectors.get_valid_state_function",
    "category": "method",
    "text": "get_valid_state_function(L, unitCellSize,, symmetries)\n\nReturns a function  states -> Bool that determines if a state is valid, i.e. satisfies the constraint imposed by the symmetry sector.\n\nArguments\n\n\'L :: Int\': length of the chain\n\'unitCellSize:: Int\' - size of the unit cell, must divide L\n\'symmetries :: Dict{String,Int}\': a dictionary of the symmetries and their sectors, e.g [\"Sz\"=> 2,\"K\"=>2]\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.make_Inversion_function-Tuple{Int64,Int64,Int64}",
    "page": "Internals",
    "title": "ED_sectors.make_Inversion_function",
    "category": "method",
    "text": "make_Inversion_function(L,Inv,a)\n\nReturns a function which computes the inversion (flip in the x-direction) of a state #Arguments\n\n\'L :: Int\': the number of sites\n\'Inv :: Int\': the Inversion sector, in {-1,1}\n\'a :: Int\': the size of the unit cell, must divide L\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.make_Z2A_function-Tuple{Int64,Int64}",
    "page": "Internals",
    "title": "ED_sectors.make_Z2A_function",
    "category": "method",
    "text": "make_Z2A_function(L,Z2A)\n\nReturns a function which computes a version of a state flipped on the B sublattice. This only makes sense when unitCellSize = 2, of course. #Arguments\n\n\'L :: Int\': the number of sites\n\'Z2A :: Int\': the Z2A sector, in {-1,1}\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.make_Z2B_function-Tuple{Int64,Int64}",
    "page": "Internals",
    "title": "ED_sectors.make_Z2B_function",
    "category": "method",
    "text": "make_Z2B_function(L,Z2B)\n\nReturns a function which computes a version of a state flipped on the B sublattice. This only makes sense when unitCellSize = 2, of course. #Arguments\n\n\'L :: Int\': the number of sites\n\'Z2B :: Int\': the Z2B sector, in {-1,1}\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.make_basis-Tuple{Int64}",
    "page": "Internals",
    "title": "ED_sectors.make_basis",
    "category": "method",
    "text": "make_basis(L; [unitCellSize=a, syms=Symmetries_Dict, constraint=constraint_function])\n\nInternal function to produce a basis from its size, constraints, and symmetry information.\n\nArguments\n\n\'L :: Int\': the number of sites\n\'unitCellSize :: Int\': the number of sites per unit cell\n\'syms :: Dict{String,Int}\': a dictionary of symmetries with their sector names, e.g. the translation sector 3 is \"K\" => 3\n\'constraint :: Function\': a function (state, L) -> bool to determine if a given state satisfies a constraint\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.make_easy_orbit_function-Tuple{Int64,Int64,Dict{String,Int64}}",
    "page": "Internals",
    "title": "ED_sectors.make_easy_orbit_function",
    "category": "method",
    "text": "make_easy_orbit_function(L,unitCellSize,symmetries)\n\nGiven the symmetry information, returns a function which computes the orbit of an individual element. Works for valid basis elements only.\n\nArguments\n\n\'L :: Int\': length of the chain\n\'unitCellSize:: Int\' - size of the unit cell, must divide L\n\'symmetries :: Dict{String,Int}\': a dictionary of the symmetries and their sectors, e.g [\"Sz\"=> 2,\"K\"=>2]\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.make_spin_flip_function-Tuple{Int64,Int64}",
    "page": "Internals",
    "title": "ED_sectors.make_spin_flip_function",
    "category": "method",
    "text": "make_spin_flip_function(L,Z2)\n\nReturns a function which computes the spin flip of a state #Arguments\n\n\'L :: Int\': the number of sites\n\'Z2 :: Int\': the Z2 sector, in {-1,1}\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.make_translation_function-Tuple{Int64,Int64,Int64}",
    "page": "Internals",
    "title": "ED_sectors.make_translation_function",
    "category": "method",
    "text": "make_translation_function(L,a,K)\n\nReturns a function which, given a state x and phase factor pf, gives the translated state T.x and new phase factor.\n\n#Arguments\n\n\'L :: Int\': the number of sites\n\'a :: Int\': the number of sites per unit cell, must divide L\n\'K :: Int\': the translation symmetry sector, 0 <= K <= L\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.measure_N_up-Tuple{UInt64,Int64}",
    "page": "Internals",
    "title": "ED_sectors.measure_N_up",
    "category": "method",
    "text": "measure_N_up(s, L)\n\nMeasures the number of spin up\'s for a state \'s\' for \'L\' spins. Returns the Sz sector of s, i.e. the number of spin\'s up.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.measure_N_up_A-Tuple{UInt64,Int64}",
    "page": "Internals",
    "title": "ED_sectors.measure_N_up_A",
    "category": "method",
    "text": "measure_N_up_A(s, L)\n\nMeasures the number of spin up\'s for a state \'s\' for \'L\' spins in the A sector, for an ABABAB sublattice structure. Returns the SzA sector of s.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.measure_N_up_B-Tuple{UInt64,Int64}",
    "page": "Internals",
    "title": "ED_sectors.measure_N_up_B",
    "category": "method",
    "text": "measure_N_up_B(s, L)\n\nMeasures the number of spin up\'s for a state \'s\' for \'L\' spins in the B sector, for an ABABAB sublattice structure. Returns the SzB sector of s.\n\n\n\n\n\n"
},

{
    "location": "internals.html#basis.jl-1",
    "page": "Internals",
    "title": "basis.jl",
    "category": "section",
    "text": "Modules = [ED_sectors]\nPages = [\"basis.jl\"]\nPublic = false\nOrder   = [:constant, :type, :function]"
},

{
    "location": "internals.html#ED_sectors.VEC",
    "page": "Internals",
    "title": "ED_sectors.VEC",
    "category": "type",
    "text": "VEC(factor,state)\n\nA type for storing a prefactor * basis_state.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.apply_S_minus-Tuple{ED_sectors.VEC,Int64}",
    "page": "Internals",
    "title": "ED_sectors.apply_S_minus",
    "category": "method",
    "text": "apply_S_minus(v,i)\n\n#Arguments *\'v: VEC\': a basis vector VEC(prefactor, state) *\'i: Int\': the index of the spin to apply the operator to, 0 <= i <= L\n\nGiven a basis vector |v>, return S_i^- |v>, which is another basis vector.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.apply_S_plus-Tuple{ED_sectors.VEC,Int64}",
    "page": "Internals",
    "title": "ED_sectors.apply_S_plus",
    "category": "method",
    "text": "apply_S_plus(v,i)\n\nGiven a basis vector |v>, return S_i^+ |v>, which is another basis vector.\n\nArguments\n\n\'v: VEC\': a basis vector VEC(prefactor, state)\n\'i: Int\': the index of the spin to apply the operator to, 0 <= i <= L\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.apply_S_x-Tuple{ED_sectors.VEC,Int64}",
    "page": "Internals",
    "title": "ED_sectors.apply_S_x",
    "category": "method",
    "text": "apply_S_x(v,i)\n\nGiven a basis vector |v>, return S_i^x |v>, which is another basis vector.\n\nArguments\n\n\'v: VEC\': a basis vector VEC(prefactor, state)\n\'i: Int\': the index of the spin to apply the operator to, 0 <= i <= L\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.apply_S_y-Tuple{ED_sectors.VEC,Int64}",
    "page": "Internals",
    "title": "ED_sectors.apply_S_y",
    "category": "method",
    "text": "apply_S_y(v,i)\n\nGiven a basis vector |v>, return S_i^y |v>, which is another basis vector.\n\nArguments\n\n\'v: VEC\': a basis vector VEC(prefactor, state)\n\'i: Int\': the index of the spin to apply the operator to, 0 <= i <= L\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.apply_S_z-Tuple{ED_sectors.VEC,Int64}",
    "page": "Internals",
    "title": "ED_sectors.apply_S_z",
    "category": "method",
    "text": "apply_S_z(v,i)\n\nArguments\n\n\'v: VEC\': a basis vector VEC(prefactor, state)\n\'i: Int\': the index of the spin to apply the operator to, 0 <= i <= L\n\nGiven a basis vector |v>, return S_i^z |v>, which is another basis vector.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.apply_operators-Tuple{ED_sectors.VEC,ED_sectors.TERM_INTERNAL}",
    "page": "Internals",
    "title": "ED_sectors.apply_operators",
    "category": "method",
    "text": "apply_operators(v,t)\n\nGiven a basis vector |v> and a term (i.e. operator) t, returns t | v>. For spin-1/2 with the X,Y,Z,+,- operators, the return vector is always another basis vector, rather than superposition thereof, which simplifies things.\n\nArguments\n\n\'v::VEC\': a vector VEC(prefactor,basis vector) to apply the term to\n\'t::TERM_INTERNAL\': a single term from a Hamiltonian to apply\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.construct_matrix_fermion_biliners-Tuple{BASIS,ABSTRACT_OP}",
    "page": "Internals",
    "title": "ED_sectors.construct_matrix_fermion_biliners",
    "category": "method",
    "text": "construct_matrix_fermion_biliners(basis, abstract_op)\n\nMakes a Hamiltonian from fermion bilinears. Only work for number-conserving Hamiltonians\n\nArgument s\n\n\'basis :: Basis\': the basis for the spin chain\n\'abstract_op :: HAMILTONIAN\': the abstract operator to implement\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.construct_matrix_full-Tuple{BASIS,ABSTRACT_OP}",
    "page": "Internals",
    "title": "ED_sectors.construct_matrix_full",
    "category": "method",
    "text": "construct_matrix_full(basis, abstract_op)\n\nAn internal method to quickly make an operator for the full basis with no symmetry constraints. \n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.construct_matrix_sym-Tuple{BASIS,ABSTRACT_OP}",
    "page": "Internals",
    "title": "ED_sectors.construct_matrix_sym",
    "category": "method",
    "text": "construct_matrix_sym(basis, abstract_op)\n\nInternal function to construct a Hamiltonian with symmetries.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.get_pauli_matrix-Tuple{String}",
    "page": "Internals",
    "title": "ED_sectors.get_pauli_matrix",
    "category": "method",
    "text": "Returns the correct Pauli matrix for the name. Really this should a hardcoded dictionary, but its essentially the same and there\'s not that many of them.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.is_Null_VEC-Tuple{ED_sectors.VEC}",
    "page": "Internals",
    "title": "ED_sectors.is_Null_VEC",
    "category": "method",
    "text": "is_Null_VEC(vec)\n\nTests if \'vec\' is the zero vector to numerical precision. Mostly used for internal testing.\n\n\n\n\n\n"
},

{
    "location": "internals.html#matrix_constructors.jl-1",
    "page": "Internals",
    "title": "matrix_constructors.jl",
    "category": "section",
    "text": "Modules = [ED_sectors]\nPages = [\"matrix_constructors.jl\"]\nPublic = false\nOrder   = [:constant, :type, :function]"
},

{
    "location": "internals.html#measurement.jl-1",
    "page": "Internals",
    "title": "measurement.jl",
    "category": "section",
    "text": "Modules = [ED_sectors]\nPages = [\"measurement.jl\"]\nPublic = false\nOrder   = [:constant, :type, :function]"
},

{
    "location": "internals.html#ED_sectors.timeseries-Tuple{Array{Complex{Float64},2},EIGSYS,Array{Float64,N} where N,Array{Complex{Float64},2},Array{Complex{Float64},2}}",
    "page": "Internals",
    "title": "ED_sectors.timeseries",
    "category": "method",
    "text": "Function to compute a two-point function Tr[rho O(t)O(0)] at a list of times {t1,t2,...}. \n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.timeseries-Tuple{Array{Complex{Float64},N} where N,EIGSYS,Array{Float64,N} where N,Array{Complex{Float64},2},Array{Complex{Float64},2}}",
    "page": "Internals",
    "title": "ED_sectors.timeseries",
    "category": "method",
    "text": "Function to compute a two-point function <psi|O(t)O(0)|psi> at a list of times {t1,t2,...}. Since we\'re doing the same operator at each time, we can optimize this a bit.\n\n\n\n\n\n"
},

{
    "location": "internals.html#ED_sectors.timeseries2",
    "page": "Internals",
    "title": "ED_sectors.timeseries2",
    "category": "function",
    "text": "Function to compute a two-point function <psi|O(t)O(0)|psi> for an eigenstate |psi> at a list of times {t1,t2,...}. Since we\'re doing the same operator at each time, we can optimize this a bit.\n\n\n\n\n\n"
},

{
    "location": "internals.html#dynamics.jl-1",
    "page": "Internals",
    "title": "dynamics.jl",
    "category": "section",
    "text": "Modules = [ED_sectors]\nPages = [\"dynamics.jl\"]\nPublic = false\nOrder   = [:constant, :type, :function]"
},

{
    "location": "internals.html#ED_sectors.display_states-Tuple{Array{UInt64,N} where N}",
    "page": "Internals",
    "title": "ED_sectors.display_states",
    "category": "method",
    "text": "Pretty printing for spin-1/2 states. #Arguments\n\n\'psi :: Array{UInt64}\': states to display\n\n\n\n\n\n"
},

{
    "location": "internals.html#utilities.jl-1",
    "page": "Internals",
    "title": "utilities.jl",
    "category": "section",
    "text": "Modules = [ED_sectors]\nPages = [\"utilities.jl\"]\nPublic = false\nOrder   = [:constant, :type, :function]"
},

]}
