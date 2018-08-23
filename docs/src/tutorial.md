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

## Overview, or Big Matrices are Hard

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

## Abstract Operators

Let's now dispense with the motivation and get down to the business of explaining how this packages works and what its important commands are.



## Choosing Bases

## Making Matrices

## Measurements 

## Dynamics