# Tutorial

```@contents
Pages = ["tutorial.md"]
Depth = 2
```

```@meta
DocTestSetup = quote
	using ED_sectors
end
```

## Overview or: Big Matrices are a Big Problem

This tutorial will introduce the basic functionality of the `ED_sectors` package. It is self-contained from a programming perspective, but assumes some knowledge of quantum mechanics.

If you're already familiar with spin chains and some of their computational difficulties, you may wish to skip to the Abstract Operators section below.

Recall that the time-independent Schrodinger Equation
```math
	\left. H |\psi_n\right> = E_n \left. |\psi_n\right>
```
determines the main physical quantities --- the energy eigenvalues and energy eigenstates of our system. Our main task is to solve this eigenvalue problem. 

For a spin-1/2 chain, the Hilbert space of ``L`` is ``\mathcal{H} = \otimes_{i=0}^{L-1} \mathbb{C}^2``, whose dimension is ``N = 2^L``. Solving the eigensystem for ``H`` naively requires ``O(N^3)`` time and ``O(N^2)`` space. In practical terms, one cannot solve this beyond about ``L = 16`` on a desktop and perhaps ``L = 26`` on the largest supercomputers, mostly due to the memory constraints. To understand a spin chain, we often want to be as close to the thermodynamic limit, ``L \to \infty`` as possible. Bigger, in this case, is better. However, the infamous exponential increase in the size of the Hilbert space limits us to very small systems. That is, unless we do something more clever.

## How to be clever

There are several ways we can be more clever and reach larger system sizes. 

### Few-Eigenstate Approaches

The simplest thing to do is to cut down on what information we need. Oftentimes, the ground state of a system is the only interest part and the rest is irrelevant; computing only a single eigenvalue is much easier than finding ``2^L``.

One can show that, starting with any state ``\left.|\phi\right>``, then repeatedly acting on it with ``e^{-\epsilon H} \approx 1 - \epsilon H`` for some small ``\epsilon > 0``, we will eventually converge to the ground state:
```math
\lim_{N \to \infty} \left( 1 - \epsilon H \right)^N \left.|\phi\right> = \left.|0\right>.
```
Physically, multiplication by ``e^{-\epsilon H}`` performs imaginary time-evolution, which cools down the system until it is stuck in its ground state. Notice, however, that all we need here is a single vector ``\left.|\phi\right>`` and the matrices we multiply it by are generically very sparse, so they are small to store and quick to multiply. This suggests that we can do much better than the naive ``O(N^3)`` expectation above.

The Lanczos Algorithm <https://en.wikipedia.org/wiki/Lanczos_algorithm> does exactly this. It gives a very quick way to find the smallest ``k`` eigenvalues and eigenvector while employing only matrix-vector multiplication rather than matrix-matrix multiplication and is therefore much, much faster than finding the full spectrum. (Actual runtime estimates are hard to find, but wikipedia estimates it as ``O(Nkd)`` where ``N`` is the size of the matrix, ``k`` eigenvectors are computed, and ``d`` is the density of the matrix.) In practice, on can find the ground states of chains of length ``L = 22`` in a few minutes on a modern desktop.

The limit here is usually storing the Hamiltonian matrix rather than the time to perform the computation. Some codes therefore implement "on-the-fly" construction of the matrix where the full matrix is never stored, but recomputed one row at a time during each matrix-vector multiplication. This can further extend the number of spins considered.

### Symmetries 

Unfortunately, there are two common situations where this is insufficient. If one cares about dynamical questions, then we often need to multiply by the propagator ``U(t) = e^{-i H t}``. If we know all the eigenvalues ``E_n`` and the (unitary) matrix of eigenvectors ``V``, then
```math
U(t) = V D(t) V^\dagger
```
where ``D(t)`` is the diagonal matrix whose elements are ``e^{-i E_n t}``. Since multiplying by diagonal matrices is quite fast, this means time-evolution is "only" a single matrix-matrix multiplication and, moreover, this is probably about the best one can do. The problem is, finding the full eigensystem is still a very hard operation, so we want to think about how to reduce its complexity as much as possible. Perhaps the best way to do this is to take advantage of symmetries of the problem.

!!! note
    An alternative approach is to give up on exact answers and instead perform approximate time-evolution. It turns out there is a clever way to do this using the Lanczos algorithm as well. Using this, one can reach times around ``t \sim 10^{4}/|J|`` where ``J`` is the energy scale of a single spin.

    This is implemented in the `quantum dynamite` package at <https://dynamite.readthedocs.io/en/latest/>. This method is particularly amenable to parallelization, enabling systems of size ``L \sim 30``.


To understand the role of symmetries, let's consider an example which will follow us through the rest of the tutorial: the XXZ spin chain:
```math
H = -\sum_{i=1}^L S_i^x S_{i+1}^x + S_i^y S_{i+1}^y + \Delta S_i^z S_{i+1}^z.
```
For the uninitiated, the XXZ chain is a model for quantum magnetism in one dimension and is one of the most-studied models in condensed matter physics. One of its crucial features is that the total spin
```math
S_{\mathrm{tot}} = \sum_i S_i^z
```
is a conserved quantity. (One can check explicitly that ``[S_\mathrm{tot}, H] = 0``.) 

This symmetry partitions out Hilbert space. The ``2^L`` states each have a representation in the ``S^z`` basis with a definite number of spin ups. For instance, the state ``\left.|\uparrow\uparrow\downarrow\uparrow\downarrow\right>`` has ``S^z = 3``. In the jargon, ``S^z`` is a "good quantum number" for this model, and has a range ``-L \le S^z \le L``. There is exactly one state (all up) with ``S^z = L``, ``L`` states with ``S^z = L-1`` and, generically, ``\binom{L}{k}`` states with ``S^z = k``, for a total of ``2^L``.

We can then think of our Hilbert space ``\mathcal{H}`` as made up of the direct *sum* of smaller Hilbert spaces ``\mathcal{H}_k`` of size ``\binom{L}{k}``, i.e.
```math
\mathcal{H} = \mathcal{H}_{-L} \oplus \cdots \mathcal{H}_{L-1} \oplus \mathcal{H}_L.
```
The fact that ``S^z`` is conserved means that the Hamiltonian act on a vector in one ``S^z``-sector ``\mathcal{H}_a`` and end up with a vector with any components in some other ``S^z``-sector ``\mathcal{H}_b``. Let's think about what this means in terms of explicit matrices. If we choose a basis for ``\mathcal{H}`` which is made by concatenating the bases for each of the ``\mathcal{H}_k`` in order, then the Hamiltonian will be block diagonal:
```math
\begin{pmatrix}
H_{-L} 	& 0 		& \cdots & 0\\
0		& H_{L-1}	& \cdots & 0\\
\vdots 	& \vdots	& \ddots & \vdots \\
0 		& 0 	& \cdots & 0 & H_{L}\\
\end{pmatrix}
```
where each block matrix ``H_k`` is the projection of ``H`` to the sector with ``k`` spins up.

In this form, it is easy to see the benefit: solving the eigenvalue problem for ``H`` is now reduced to solving the problem for each of the ``H_k``'s. For ``L = 24``, the size of the largest block is ``\binom{24}{12} \approx 2^{21}``, a substantial reduction in the difficulty of the problem!

Here we worked with the particular example of the ``S_{\mathrm{tot}}`` symmetry, but the idea is general. Every (Abelian) symmetry of the problem can be used to further reduce the sizes of the diagonal blocks involved. The most common other symmetries are translation, parity, and inversion. Using some combination of these, one can often reduce the effective size of the system by 2 to 5 spins, and sometimes more. I believe the world record with this method is around ``L = 48``. 

While a few extra spins might not seem like a lot, consider that most dynamical simulations max out at ``L = 14``. In comparison, working with ``L = 20`` is a huge improvement!

The downside of this method is that actually putting the matrix into a block-diagonal form requires a lot of bookkeeping. In fact, the ``S_{\mathrm{tot}}`` example was chosen above since it is particularly simple. Generically, finding the block-diagonal form of the matrix involves a fairly complicated and non-local change of basis. Its form is dictated by the (projective) representations of the finite symmetry groups. 

!!! note
    **The purpose of this package is to abstract away this troublesome procedure. With `ED_sectors`, finding a single symmetry block of the Hamiltonian can be done in a single line.**

## Defining Hamiltonians: Abstract Operators

Let's now dispense with the motivation and get down to the business of explaining how this packages works and what its important commands are. We will start with defining Hamiltonians. For a spin chain of length ``L``, the space of Hamiltonians has dimension ``4^L`` -- its huge! Physically interesting Hamiltonians, however, are only a small corner of this space. They tend to be translation invariant, or at least nearly so, and involve only local or quasi-local interactions. In practice, then, they are usually a sum of a few terms. Similarly, most observables involve only a few sites, or perhaps a sum of the same few-site operator translated across the system. In `ED_sectors` we adopt a syntax designed to make writing physical Hamiltonians and common observables as easy as possible.

To do this, we have two main data structures: [`ABSTRACT_OP`](@ref)'s and [`TERM`](@ref)'s. The former, short for "abstract operator", is a representation of an operator as a string of Pauli matrices and, just how a Hamiltonian is written as a sum of individual terms, an `ABSTRACT_OP` is made by adding `TERM`s together.

As an example, let's define the XXZ Hamiltonian mentioned above:
```julia
XXZ = ABSTRACT_OP(8)
XXZ += TERM("XX")
XXZ += TERM("YY")
XXZ += 0.3*TERM("ZZ")
```
In the first line, we create an abstract operator with `16` sites and, by default, periodic boundary conditions. We then add three terms to the Hamiltonian which are automatically translated across all sites.

`ABSTRACT_OP`s have a choice of either periodic or open boundary conditions, which is easily adjusted with an optional argument. Additionally, we can give abstract operators names to make them easier to debug. For detailed documentation, see [`ABSTRACT_OP`](@ref).

```
ising = ABSTRACT_OP(15; name="Ising Model", pbc=false) + TERM("ZZ") + TERM("X")
```

Just like real Hamiltonians, we can add abstract operators together, or multiply them by a (complex) prefactor. For instance, we can easily define a magnetic field in the ``Z`` direction, adjust its magnitude, and add it on to our XXZ operator.
```
magnetic_field = ABSTRACT_OP(8) + TERM("Z")
magnetic_field_2 = 0.1*magnetic_field
XXZ_B = XXZ + magnetic_field_2
```

This syntax is great for creating relatively complex operators, like Hamiltonians. But individual observables are often much simpler and often live only on a single site or a few sites. For this situation, there is a notational shortcut.
```julia
middle_field = ABSTRACT_OP(8,"Z",4)
```
This creates the operator ``S_4^z`` in a single line. Alternatively, a single `TERM` can be supplied instead. With this we can easily define, for instance, an order parameter ``\tfrac{1}{L} \sum_i S_i^z``.
```julia
magnetic_order_parameter = ABSTRACT_OP(8,(1/8)*TERM("Z"))
```



To create all the different types of operators you might encounter in the course of your day-to-day spin chain exploration, there are several different ways to make `TERM`'s. Some examples are provided here, and all the options are spelled out at [`TERM`](@ref).

| Example                                     | Result                                                                                           |
|:-------------------------------------------:|:------------------------------------------------------------------------------------------------:|
| `TERM("ZZ")`                                | ``\displaystyle \sum_{i=0}^{L-1} S_i^z S_{i+1}^z``                                               |
| `TERM("ZIZ",4)`                             | ``S_4^z S^z_6 ``                                                                                 |
| `TERM(3.0,"X",3)`                           | ``\displaystyle3 \sum_{i=3}^{L-1} S_i^x ``                                                       |
| `TERM("YY",3,repeat=2)`			          | ``\displaystyle\sum_{\substack{i=2\\ i \mathrm{ even}}} S_i^y S_{i+1}^y``                        |
| `TERM(4.3,"ZXZ",5,repeat=3)`                | ``\displaystyle 4.3 \sum_{\substack{i=5\\ i \equiv 0 \pmod 3}}^{L-1} S_i^z S_{i+1}^x S_{i+2}^z ``|
| `TERM(0.5,Dict(1=>"X", 7=>"X))`             | ``S_1^x S_7^x``                                                                                  |
| `TERM(0.5,Dict(1=>"X", 7=>"X); repeat = 1)` | ``\displaystyle \sum_{i=1}^{L-1} S_i^x S_{i+7}^x``                                               |

The several different ways of making `TERM`'s shown here are general enough to make any Hamiltonian, even if it has next-next-nearest neighbor interactions, isn't translationally invariant, or other unusual possibilities. If the most convenient syntax is insufficient then, as shown in the last two examples, one can define terms based on a dictionary `Dict(site => "op")`. Adding these together can make any operator whatsoever.


!!! warning "Index Convention"
    Abstract operators are zero-indexed. The left-most spin is spin ``0`` and, in a chain of length ``L``, the right-most spin has index ``L-1``.

!!! warning "Boundary Conditions"
    The above examples use periodic boundary conditions. In that case, all indices larger than ``L-1`` will wrap around, i.e. indices should be interpreted ``\pmod L``.

    With open boundary conditions, however, all terms whose indices exceed ``L-1`` will simply be dropped. 


## Symmetries and Bases

Making a [`BASIS`](@ref) is as easy as specifying the number of sites. Optionally, one can also include the size of the unit cell, symmetries, and a constraint on which sites are allowed in the Hilbert space.


For instance, making a basis for `16` spins is just
```julia
BASIS(16)
```

As discussed extensively above, however, the real purpose of making a basis in exact diagonalization is to take symmetries into account. There are a number of symmetries built-in, which we now list in a convenient table.


| Name       | Symmetry                             | Values                  | Generator                                                                           |
|:----------:|--------------------------------------|:-----------------------:| ------------------------------------------------------------------------------------|
| Tr         | Translation                          | ``1 \le t \le L/a``     | ``S_i^\mu \to S_{i+t \pmod{L/a}}``                                                   |
| P          | Parity (Spin flip)                   | ``\{-1,1\}``            | ``\displaystyle \prod_{i=0}^{L-1} S_i^x``                                           | 
| I          | Inversion                            | ``\{-1,1\}``            | ``S_i^\mu \to S_{L-1-i}^\mu``                                                       |
| Sz         | Total ``S^z``                        | ``-L \le Sz \le L``     | n/a; counts ``\sum_{i=0}^{L-1} S_i^z``                                              |
| PA         | P for even spins                     | ``\{-1,1\}``            | ``\displaystyle \prod_{\substack{i=0\\ i \equiv 0 \pmod 2}}^{L-1} S_i^x``           |
| PB         | P for odd spins                      | ``\{-1,1\}``            | ``\displaystyle \prod_{\substack{i=0\\ i \equiv 1 \pmod 2}}^{L-1} S_i^x``           |
| SzA        | ``S^z`` for even spins               | ``-L/2 \le SzA \le L/2``| n/a; counts ``\displaystyle\sum_{\substack{i=0\\ i \equiv 0 \pmod 2}}^{L-1} S_i^z`` |
| SzB        | ``S^z`` for odd spins                | ``-L/2 \le SzB \le L/2``| n/a; counts ``\displaystyle\sum_{\substack{i=0\\ i \equiv 1 \pmod 2}}^{L-1} S_i^z`` |

To make a basis for the third translation sector and odd parity sector for instance, is simply
```julia
basis_sym = BASIS(16; syms=Dict("Tr" => 3, "P" => -1))
```

One can also add custom constraints onto the Hilbert space. Suppose, for instance, you wanted to use the "Rydberg constraint" that two adjacent spins cannot both point up at once. After designing a function to test for that, it can be easily used to make bases satisfying those constraints.
```julia
function test_Rydberg(state) :: Bool
	[...]
end

basis_Rydberg = BASIS(24; constraint = test_Rydberg)
```

More information about making bases can be found at `BASIS`](@ref).

## Making Actual Matrices

Making actual matrices is done by specifiying the abstract operators, and the basis where it should live.

```julia
L = 16
ising = ABSTRACT_OP(L) + TERM("ZZ") + 0.3*TERM("X")
basis_sym = BASIS(L; syms=Dict("Tr" => 3, "P" => -1))

H = Operator(ising, basis_sym)
```
This produce a sparse matrix `H` which we can diagonalize or otherwise manipulate as we like.

There is also a shortcut syntax for making single-site or single-term observables.
```julia
Sx2 = Operator("X",2,basis)
```


!!! warning "Invalid Symmetries"
    `ED_sectors` does not (yet!) check if an operator actually *has* a particular symmetry. If the operator doesn't actually have the symmetry you specify, you will get nonsense!


## Measurements 

Now that we are equipped with the tools for making matrices, we can compute expectation values to perform measurements. At this point, this is basically just linear algebra, but several convenient wrapper functions are provided to speed things up.

As an example, let's measure the magnetization on the third spin in the ground state 
```julia
L = 16
ising = ABSTRACT_OP(L) + TERM("ZZ") + 2.3*TERM("X")
basis_sym = BASIS(L; syms=Dict("Sz" => 0))

H = Operator(ising, basis_sym)

(E_0, psi_0) = ground_state(H)

Sz3 = Operator("Z",3,basis_sym)

expectation(Sz3,psi_0)
```

There are also more complex options, such as thermal density matrices. Computing them, of course, requires finding the full spectrum. This is done with the [`EIGSYS`](@ref) command and, as one might expect, is quite computationally expensive. 

```julia
L = 16
ising = ABSTRACT_OP(L) + TERM("ZZ") + 2.3*TERM("X")
basis_sym = BASIS(L; syms=Dict("Sz" => 0))

H = Operator(ising, basis_sym)
Sz3 = Operator("Z",3,basis_sym)

full_spectrum = EIGSYS(H)

beta = 3.0
rho_beta = thermal_density_matrix(full_spectrum,beta)

expectation(Sz3,rho_beta)
```

These examples are just a small taste of the measurements one can easily perform. For more information on measurements, see the specific examples in the section on [Ground State Measurements](measurement.md), [Full Spectrum Measurements](full_spec.md), or [Dynamics](dynamics.md). In particular, there are easy functions to compute entanglement entropy and time-evolution of both states and operators.

## Conclusion

Thanks for reading the tutorial! As you probably realized, this guide is a little front-heavy. This is because most of the measurement and dynamics functions are almost exactly the same as their mathematical counterparts, so they should be relatively self-explanatory from their documentation.

This package is still a work-in-progress, so there are several unfinished or undocumented features hiding in the code. Here's the roadmap for things that will be added soon:

- Other types of sites, like spinless fermions or spin-1 particles
- Custom or user-defined symmetries
- Checking that an operator satisfies the given symmetries for a basis
- Autocorrelations at temperature with more efficient code.
- Linking to dynamite?
- Projectors from symmetry sectors into the full Hilbert space, which are useful for computing entanglement and doing dynamics with operators that don't respect the symmetries.

For any comments, questions, or issues, please start a GitHub issue.

