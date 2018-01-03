    #!/usr/bin/python2.7

from time import time
import math as m, numpy as np

import scipy.sparse.linalg as la
from scipy.sparse import *
from scipy import io 

import scipy.optimize as opt

import numpy.linalg as lafull 




#--------------------construct Hamiltonians, etc


pauli_1 = csr_matrix(np.array([[1, 0], [0, 1]]))
pauli_x = csr_matrix(np.array([[0, 1], [1, 0]]))
ipauli_y = csr_matrix(np.array([[0, 1], [-1, 0]])) #keeping it real
pauli_z = csr_matrix(np.array([[1,0], [0, -1]]))


xx = kron(pauli_x,pauli_x).asformat("csr")
yy = ((-1.0) * kron(ipauli_y,ipauli_y)).asformat("csr")
zz = kron(pauli_z,pauli_z).asformat("csr")

x1x = kron(pauli_x,kron(pauli_1,pauli_x)).asformat("csr")
x1x.eliminate_zeros()

y1y = (-1.0)*kron(ipauli_y,kron(pauli_1,ipauli_y)).asformat("csr")
y1y.eliminate_zeros()

z1z = kron(pauli_z,kron(pauli_1,pauli_z)).asformat("csr")
z1z.eliminate_zeros()

zxz = kron(pauli_z,kron(pauli_x,pauli_z)).asformat("csr")
zxz.eliminate_zeros()

zyyz = -kron(pauli_z,kron(ipauli_y,kron(ipauli_y,pauli_z))).asformat("csr")
zyyz.eliminate_zeros()

zzzz = kron(pauli_z,kron(pauli_z,kron(pauli_z,pauli_z))).asformat("csr")
zzzz.eliminate_zeros()

zx1xz = kron(pauli_z,kron(pauli_x,kron(pauli_1,kron(pauli_x,pauli_z)))).asformat("csr")
zx1xz.eliminate_zeros()

zy1yz = (-1.0)*kron(pauli_z,kron(ipauli_y,kron(pauli_1,kron(ipauli_y,pauli_z)))).asformat("csr")
zy1yz.eliminate_zeros()

#n-particle interaction embedding into the whole Hilbert space
def shift(mat,left,right):
    L = eye(2**left)
    R = eye(2**right)
    M = kron(kron(L,mat),R)
    return M

def sigma_X_sym(L):
    U = eye(2**L)
    for i in range(0,L-1,2):
        flipSpin = shift(pauli_x,i,L-i-1)
        U = U * flipSpin
    return U

def tau_X_sym(L):
    U = eye(2**L)
    for i in range(1,L,2):
        flipSpin = shift(pauli_x,i,L-i-1)
        U = U * flipSpin
    return U

def X_sym(L):
    U = eye(2**L)
    for i in range(1,L):
        flipSpin = shift(pauli_x,i,L-i-1)
        U = U * flipSpin
    return U

def make_H_Ising(L,g_sigma,gamma):
    H = csr_matrix((2**L,2**L))
    for i in range(0,L-2,2):
        H = H - shift(zz,i,L-i-2) #Z Z
        H = H - shift(pauli_x,i,L-i-1)    #X
        #H = H - gamma*shift(xx,i,L-i-2) #X_A X_B


    wrap_mat_ZZ = kron(pauli_z,kron(eye(2**(L-2)),pauli_z)).asformat("csr")
    wrap_mat_ZZ.eliminate_zeros()
    #wrap_mat_XX = kron(pauli_x,kron(eye(2**(L-2)),pauli_x)).asformat("csr")
    #wrap_mat_XX.eliminate_zeros()
    H = H - wrap_mat_ZZ #- gamma*wrap_mat_ZZ

    return H

#make the Hamiltonian matrix with periodic boundary conditions
def make_H_g_Trivial(L,g_sigma,gamma,beta,mag):
    #L          - Length
    #J          - A ZZ coupling
    #h          - A X strength
    #B_scale    - overall scale for the B lattice
    #g          - B ZZ coupling 'interaction'
    #mag        - magnetic field strength
    #mag_pos    - site for the magnetic field from 0 to (L-1)
    # B_scale = 10
    # g_tau = 0.2
    # #g_sigma = 1
    # U_tau = 0.1
    # U_sigma = 0.2
    B_scale = 1
    g_tau = 0.0
    #g_sigma = 1
    U_tau = 0.0
    U_sigma = 0.0

    #mag = 0.001

    H = csr_matrix((2**L,2**L))
    for i in range(0,L-2,2):
        H = H - g_sigma*shift(z1z,i,L-i-3) #Z_A 1 Z_A
        H = H - U_sigma*shift(x1x,i,L-i-3) #X_A 1 X_A
    wrap_mat_Z1Z = kron(pauli_z,kron(eye(2**(L-3)),pauli_z)).asformat("csr")
    wrap_mat_Z1Z.eliminate_zeros()
    H = H - g_sigma*shift(wrap_mat_Z1Z,0,1)

    wrap_mat_X1X = kron(pauli_x,kron(eye(2**(L-3)),pauli_x)).asformat("csr")
    wrap_mat_X1X.eliminate_zeros()
    H = H - U_sigma*shift(wrap_mat_X1X,0,1)

    for i in range(0,L,2):
        #print "X_A " + str(i)
        H = H - shift(pauli_x,i,L-i-1)

    for i in range(1,L,2):
        #print "X_B " + str(i)
        H = H -  B_scale*shift(pauli_x,i,L-i-1) #X_B

    for i in range(1,L-2,2):
        H = H - B_scale*g_tau*shift(z1z,i,L-i-3) #Z_B 1 Z_B
        H = H - B_scale*U_tau*shift(x1x,i,L-i-3) #Z_B 1 Z_B
    # H = H - B_scale*g_tau*shift(wrap_mat_Z1Z,1,0)
    # H = H - B_scale*U_tau*shift(wrap_mat_X1X,1,0)

    for i in range(0,L-1,2):
        #print "X_AX_B " + str(i) + " to " + str(i+2-1)
        H = H - gamma*shift(xx,i,L-i-2) #X_A X_B

    for i in range(0,L-4,2):
        #print "ZZZZ " + str(i) + " to " + str(i+4-1)
        H = H - beta*shift(zzzz,i,L-i-4) #X_A X_B    

    for i in range(0,L):
        H = H - mag*shift(pauli_z,i,L-i-1)

    return H


def make_H_g_SPT(L,g_sigma,gamma,Bgap,mag):
    #adjusts the B gap
    B_scale = Bgap
    #adds X terms to the first and last A and B on each side
    x_edge = 0.2
    #a small magnetic field in the Z direction
    #mag = 0.001

    #the ZZ coupling for the A spins
    #g_sigma = 1
    #the ZZ coupling for the B spins, to move away from the zero-correlation length paramagnet
    g_tau = 0.2
    #XX interation for the A spins
    U_sigma = 0.2
    #XX interaction for the B spins
    U_tau = 0.1
    #X_A X_B interaction --- obeys symmetries
    #gamma = 0#0.1
    #Z_AZ_B Z_A Z_B interaction --- obeys symmetries
    beta = 0#0.15

    H = csr_matrix((2**L,2**L))
    for i in range(0,L-2,2):
        H = H - g_sigma*shift(z1z,i,L-i-3) #Z_A 1 Z_A

    for i in range(0,L-2,2):
        H = H - B_scale*shift(zxz,i,L-i-3) #X_B

    for i in range(1,L-2,2):
       H = H - shift(zxz,i,L-i-3)#X_A
       H = H - B_scale*g_tau*shift(z1z,i,L-i-3) #Z_B 1 Z_B

    #interaction for Ising
    for i in range(1,L-4,2):
        H = H - U_sigma*shift(zx1xz,i,L-i-5) #X_A 1 X_A

    #interaction for B's
    for i in range(0,L-4,2):
        print("ZX1XZ " + str(i) + " to " + str(i+5-1))
        H = H - B_scale*U_tau*shift(zx1xz,i,L-i-5) #X_B 1 X_B

    #edge fields to break the edge qbits
    H = H - x_edge*shift(pauli_x,0,L-1) - x_edge*B_scale*shift(pauli_x,1,L-2)
    H = H - x_edge*B_scale*shift(pauli_x,L-1,0) - x_edge*shift(pauli_x,L-2,1) 

    #H = H - 0.2*shift(x1x,0,L-3)


    for i in range(1,L-3,2):
        #print "ZYYZ " + str(i) + " to " + str(i+4-1)
        H = H - gamma*shift(zyyz,i,L-i-4) #X_A X_B

    for i in range(0,L-4,2):
        print("ZZZZ " + str(i) + " to " + str(i+4-1))
        H = H - beta*shift(zzzz,i,L-i-4) #X_A X_B    

    for i in range(0,L):
    	H = H - mag*shift(pauli_z,i,L-i-1)
    #H = H - mag*shift(pauli_z,0,L-1)

    return H


def get_lower_spectrum(H,num):
    (w,v) = la.eigsh(H, k=num, which='SA')
    return w

def get_GS(H,kk):
    (w,v) = la.eigsh(H, k=kk, which='SA')
    ind_GS = np.argmin(w)
    GS = v[:,ind_GS]
    return GS

#make magnetization matrix at pos in [0,L-1]
def mag(pos,L):
    return shift(pauli_z,pos,L-1-pos)


#make magnetization matrix at pos in [0,L-1]
def magX(pos,L):
    return shift(pauli_x,pos,L-1-pos)


#make magnetization matrix at pos in [0,L-1]
def magY(pos,L):
    return shift(pauli_y,pos,L-1-pos)


#make magnetization matrix at pos in [0,L-1]
def magZ(pos,L):
    return shift(pauli_z,pos,L-1-pos)


def getLGapDifference(L2,L1,g_sigma,gamma,beta):
    H2 = make_H_g_SPT(L2,g_sigma,gamma,beta)
    H1 = make_H_g_SPT(L1,g_sigma,gamma,beta)
    (w1,v1) = la.eigsh(H1, k=5, which='SA')
    (w2,v2) = la.eigsh(H2, k=5, which='SA')
    LE1 =  L1*(w1[1]-w1[0])
    LE2 =  L2*(w2[1]-w2[0])
    return LE2-LE1

#switch to using a built-in Newton's method, it's much faster...
def transition_finder2(L1,L2,guess,gamma,beta):
    zero = opt.newton_krylov(lambda g:getLGapDifference(L2,L1,g,gamma,beta),guess,verbose=1)
    return zero


# numEvs=20
# g_sigma = 1.421
# gamma = 0.1
# Lmin = 8
# Lmax = 18
# print "g_sigma: "+ str(g_sigma)
# for L in range(Lmin,Lmax,4):
#     print "L: "+ str(L)
#     H = make_H_g_SPT(L,g_sigma,gamma,0)
#     (w,v) = la.eigsh(H, k=numEvs, which='SA')
#     print v.shape
#     magnetization =  [[np.dot(v[:,k],mag(i,L) * v[:,k]) for i in range(0,L)] for k in range(0,20)]
#     np.savetxt("data/Ising_tw_AB_interaction_"+str(L)+"_gsigma_"+str(g_sigma)+".csv",w,delimiter=",")
#     np.savetxt("data/Ising_tw_AB_interaction_mag_"+str(L)+"_gsigma_"+str(g_sigma)+".csv",magnetization,delimiter=",")


def getExpForDiagonal(spec,t):
    expspecs = [np.exp(1j*t*s) for s in spec]
    print(t)
    return np.diag(expspecs)


#symmetry matrix
def G_momentum(Ly):
    G = dok_matrix((2**Ly,2**Ly), dtype='i')
    for k in range(0,2**(Ly-1)):
        G[k,2*k] = 1
        G[k+2**(Ly-1),1+2*k] = 1
    return G.asformat("csr")



#quick way to produce the reduced density matrix
#l is the length of the partition we want to have left
#L is the total number of sites
def reduce_density_matrix(phi, l, L):
    dim = 2**l# dimension of the reduced space
    phi_split = np.split(phi,dim)
    rho_reduced = np.ndarray(shape=(dim,dim),dtype='complex128')
    for x in range(0,dim):
        for y in range(0,dim):
            rho_reduced[x,y] = np.dot(np.conjugate(phi_split[x]),phi_split[y])
    return rho_reduced



#get the reduced density matrix for the edge spins
#for a  state (which has probably been time evolved the way I'm using this)
def reduce_density_matrix_edge(phi_t,L):
    # phi2 = np.ravel(np.dot(G_momentum(L).todense(),phi_t))
    phi2 = np.ravel(np.dot(np.matmul(G_momentum(L).todense(),G_momentum(L).todense()),phi_t))
    rhoAYZ = reduce_density_matrix(phi2,3,L)
    # rhoAYZ2 = np.reshape(rhoAYZ,(2,2,2,2,2,2))
    # print(rhoAYZ2.shape)
    # rhoAY = np.ndarray(shape=(4,4))
    # # for x in range(0,4):
    # #     for y in range(0,4):
    # #         rhoAY[x,y] = 
    # print(rhoAYZ2[1,1])
    return rhoAYZ

L = 12
g_sigma = 1.0
gamma = 0
beta = 0
mag = 0
H = make_H_g_Trivial(L,g_sigma,gamma,beta,mag)
(w,v) = la.eigsh(H, k=30, which='SA')
print(w)

# L = 6

# H = make_H_Ising(L,1,0) 

# (w,v) = la.eigsh(H, k=10, which='SA')

# H2=H.todense()
# (spec, vecs) = lafull.eigh(H2)
# print("Diagonalized")

# print(w)
# print(vecs.shape)
# print(v[:,0].shape)

# g = G_momentum(L)
# x = X_sym(L)
# for k in range(0,10):
# 	a = v[:,k]
# 	b = g*a
# 	c = x*a
# 	d = H*a
# 	i = 0
# 	while np.abs(a[i]) < 10e-6:
# 		print(np.abs(a[i]))
# 		i = i+1

# 	print(i,a[i])
# 	print("H",k,d[i]/a[i])
# 	print("T",k,b[i]/a[i])
# 	print("X",k,c[i]/a[i])
# 	print()


# numEvs=20
# gamma = 0
# Lmin = 12
# Lmax = 13
# beta = 0.1
# # g_sigma= 0.1
# # g_sigma_disp = format(g_sigma,'.3f')

# #print(transition_finder2(8,12,1.0,0,0))



# #g_sigma = 1.0

# for g_sigma in [1.0]:#Bgap in np.arange(10,11,1): #[0.1,0.8,1.309,1.5,2.0]:
#     #xfield_disp = format(xfield,'.4f')
#     #g_sigma = 1.0
#     Bgap = 5
#     Bgap_disp = format(Bgap,'.4f')
#     beta_disp = format(beta,'.4f')
#     g_sigma_disp = format(g_sigma,'.4f')
#     print("g_sigma: "+ str(g_sigma_disp))
#     for L in range(Lmin,Lmax,4):
#         print("L: "+ str(L))
#         print(g_sigma)
#         H = make_H_g_SPT(L,g_sigma,0,Bgap,0.01)
#         (w,v) = lafull.eigh(H.todense())
#         #(w,v) = la.eigsh(H, k=32, which='SA')
#         print(w[0])
#         #vv = v[:,0]
#         #print(lafull.norm(H*vv-w[0]*vv))



#         magnetization =  np.array([[np.dot(np.transpose(v[:,k]),magZ(i,L) * v[:,k])[0,0] for i in range(0,L)] for k in range(0,2**L)])
#         print("mag",magnetization)
#         # np.savetxt("dataIsing/Ising_tw_AB_interaction_"+str(L)+"_gsigma_"+str(g_sigma_disp)+".csv",w,delimiter=",")
        
#         #np.savetxt("dataIsing4/mag_L_"+str(L)+"_gsigma_"+str(g_sigma_disp)+".csv",magnetization,delimiter=",")

#         # print(sigma_X_sym(8).shape)

#         # v2=np.transpose(v[:,1])
#         # print(np.dot(v2,np.matmul(magZ(0,L).todense(),v[:,0])))

#         H_dynamic = make_H_g_SPT(L,g_sigma,0,Bgap,0)
#         H2=H_dynamic.todense()
#         (spec, vecs) = lafull.eigh(H2)
#         print("Diagonalized")



#         #print(spec[:10])

#         O = np.array(vecs,copy=True)
#         Ot = np.transpose(vecs)

#         D = np.diag(spec)
#         D2 = np.diag(spec*spec)

#         Sz=magZ(0,L).todense()
#         #Sz=magZ(2,L).todense()

#         OtSz= np.matmul(Ot,Sz)
#         OtSzO = np.matmul(OtSz,O)

#         Z = np.sum([np.exp(-beta*w[n]) for n in range(len(spec))])
#         print("Z:"+str(Z))



#         ts=np.logspace(-2,6,num=20)
#         states = np.transpose(v)
#         print(len(states))
#         coherences3 = []
#         for t in ts:
#             correlation_mat = np.matmul(O,
#                     np.matmul(getExpForDiagonal(spec,-t),
#                     np.matmul(OtSzO,
#                     np.matmul(getExpForDiagonal(spec,t),OtSz))))
#             #print("time: "+str(t))
#             #print([np.dot(st,np.matmul(correlation_mat,np.transpose(st)))[0,0] for st in states])
#             corrSum = (1.0/Z)*np.sum([np.exp(-beta*w[n])*np.dot(states[n],np.matmul(correlation_mat,np.transpose(states[n])))[0,0] for n in range(0,len(states))])
#             corr2 = np.array([t,np.real(corrSum),np.imag(corrSum)])
#             coherences3.append(corr2)
#         coherences4 = np.array(coherences3)
#         #print(coherences4)



#         np.savetxt("dataIsing4/coherences_sum_L_"+str(L)+"_gsigma_"+str(g_sigma_disp)+"_beta_"+beta_disp+"_Bgap_"+Bgap_disp+".csv",coherences4,delimiter=",")

#        #  rho_ts = []
#        #  for t in ts:
#        #      phi_t = np.dot(
#        #              np.matmul(O,
#        #              np.matmul(getExpForDiagonal(spec,-t),
#        #              Ot)),
#        #              np.transpose(states[0]))
#        #      # phi_t = np.transpose(states[0])
#        #      phi_t = np.ravel(phi_t)
#        #      rdm = reduce_density_matrix_edge(phi_t,L).flatten().view(float)
#        #      #rdm = reduce_density_matrix(phi_t,1,L).flatten().view('float64')
#        #      rdm2 = np.concatenate((np.array([t]),rdm))
#        #      rho_ts.append(rdm2)
#        #  rho_ts2 = np.asarray(rho_ts)
#        # #print(rho_ts2)

#        #  np.savetxt("dataIsing4/rdm_GS_L_"+str(L)+"_gsigma_"+str(g_sigma_disp)+".csv",rho_ts2,delimiter=",")

#         #print([[(i,j) for i in range(3) for j in range(3)]])