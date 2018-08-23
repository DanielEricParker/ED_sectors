# ED_sectors Documentation

Welcome to `ED_sectors`, a package which provides a fast and straightforward interface for doing exact diagonalization. Its main purpose is to diagonalize Hamiltonians for spin-1/2 chains, taking advantage of symmetries to speed up the computations -- hence "sectors". `ED_sectors` is written in Julia, a new language that combines python-like syntax with speed comperable to C and was designed specifically for scientific computation. The code is designed to have easy-to-use syntax while remaining speedy.


## Quick Example

Let's start off with a quick example of an `ED_sectors` program just to see how easy it can be. Let's take one of the most popular spin chains, the XXZ model:
```math
H = -\sum_{i=1}^L S_i^x S_{i+1}^x + S_i^y S_{i+1}^y + \Delta S_i^z S_{i+1}^z.
```
This has several symmetries, including translation and conservation of total spin ``S_{\mathrm{tot}} = \sum_i S_i^z``. We can easily take advantage of those to reduce the dimension of the Hilbert space.

```julia
L = 16 #length of the spin chain

#create the basis with total Sz = +21and work in the second translation sector
basis = BASIS(L; syms = Dict("Sz" => 1, "Tr" => 2))

#create an abstract representation of the XXZ Hamiltonian
XXZ = ABSTRACT_OP(L; name = "XXZ", pbc = true)
XXZ -= TERM("XX") - TERM("YY") - 0.3TERM("ZZ")

#construct an explicit matrix for the operator in the chosen basis
H = Operator(XXZ, basis)

#get the ground state energy and wavefunction
(E_0, psi_0) = k_eigvals(H,1)

#compute the magnetization in the ground state

S_tot = Operator(TERM("Z"),basis)

mag_Z = expectation(S_tot,psi_0)
```

In only a few lines we can define a basis, create a Hamiltonian, find its ground state, and measure expectation values. See the Tutorial or the Full API to learn more!

## Installation

Eventually `ED_sectors` will be a true Julia package, but for now you need to install it manually. To install ED_Sectors:
1. Install Julia 1.0 or later from <https://julialang.org/downloads/>. Linux repositories often have older versions of Julia, so it is advised that you use the latest version from the Julia website. Installing Julia is painless: you just drop a folder onto your computer. No linking libraries is necessary!

2. Clone this package by navigating to a directory of your choice and running  
   ```
   git clone https://github.com/DanielEricParker/ED_sectors
   ```

3. Download and install the necessary Julia libraries by running the following commands in the Julia REPL.  
   ```julia
   using Pkg
   Pkg.add("Nullables")
   Pkg.add("IterativeSolvers")
   Pkg.add("Arpack")
   ```

4. Include `ED_sectors` in your julia program via  
   ```julia
   push!(LOAD_PATH,"[your path to ED_sectors here]/src/")
   using ED_sectors
   ```

## Contents

To learn about the "moving pieces" in `ED_sectors`, please see the Tutorial. More advanced user may wish to consult the Full API documentation. Many common computations are already built-in, including computing entanglement entropy, autocorrelations in time, measurements at temperature, and finding phase transitions! In the case you need to extend the functionality, the internal functions and datatypes are documented under Internals.


```@contents

```

```@meta
DocTestSetup = quote
	using ED_sectors
end
```


This package was created by Daniel E. Parker at UC Berkeley. Special thanks to Thomas Scaffidi and Romain Vasseur for inspiring me to create this package and testing out much of its functionality.