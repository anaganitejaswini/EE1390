import numpy as np 
import sympy as sp 
import matplotlib.pyplot as plt
import subprocess
import shlex
valOfK = 0
k = sp.symbols('k',integer=True)
solutions = []
areaValue = 28 ;

A = np.array([k,-3*k,1])
B = np.array([5,k,1])
C = np.array([-k,2,1])
	

areaMatrix = np.vstack((A,B,C))

areaMatrix = sp.Matrix(areaMatrix)
areaMatrixDeterminant1 = areaMatrix.det() - 56
areaMatrixDeterminant2 = areaMatrix.det() + 56
solutions.append(sp.solve(areaMatrixDeterminant1, k))
solutions.append(sp.solve(areaMatrixDeterminant2, k))
print(solutions)
valOfK = solutions[0][0]
print(valOfK)
k1 = int(valOfK)
A = np.array([k1,-3*k1])
B = np.array([5,k1])
C = np.array([-k1,2])
def proj(A,B,C):
	 X = ((A - B).T)
	 Y = ((C - B))
	 F = np.matmul(X,Y)
	 G = np.matmul(Y.T,Y)
	 D = ((F/G)*(C - B)) + B 
	 return D
def line_intersect(AD,BE):
	n1 = (B - C).T
	n2 = (A - C).T
	N = np.vstack((n1,n2))
	p = np.zeros(2)
	p[0] = np.matmul(n1,A)
	p[1] = np.matmul(n2,B)
	return np.matmul(np.linalg.inv(N),p)
BC = np.vstack((B,C)).T
CA = np.vstack((C,A)).T
D = proj(A,B,C)
E = proj(B,C,A)
AD = np.vstack((A,D)).T
BE = np.vstack((B,E)).T

print (line_intersect(AD,BE))
H = line_intersect(AD,BE)
len = 10
lam_1 = np.linspace(0,1,len)
x_AB = np.zeros((2,len))
x_BC = np.zeros((2,len))
x_CA = np.zeros((2,len))
x_AD = np.zeros((2,len))
x_BE = np.zeros((2,len))

for i in range (len):
	temp1 = A + lam_1[i]*(B - A)
	x_AB[:,i] = temp1.T
	temp2 = B + lam_1[i]*(C-B)
	x_BC[:,i] = temp2.T
	temp3 = C + lam_1[i]*(A-C)
	x_CA[:,i]=temp3.T
	temp4 = A + lam_1[i]*(D-A)
	x_AD[:,i]=temp4.T
	temp5 = B + lam_1[i]*(E-B)
	x_BE[:,i]=temp5.T

plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')

plt.plot(A[0],A[1],'o')
plt.text(A[0]*(1+0.1),A[1]*(1-0.1),'A')
plt.plot(B[0],B[1],'o')
plt.text(B[0]*(1-0.2),B[1]*(1),'B')
plt.plot(C[0],C[1],'o')
plt.text(C[0]*(1+0.03),C[1]*(1-0.1),'C')
plt.plot(D[0],D[1],'o')
plt.text(D[0]*(1+0.1),D[1]*(1-0.1),'D')
plt.plot(E[0],E[1],'o')
plt.text(E[0]*(1+0.1),E[1]*(1-0.1),'E')
plt.plot(H[0],H[1],'o')
plt.text(H[0]*(1+0.1),H[1]*(1-0.1),'H')


plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc = 'best')
plt.grid()
plt.savefig('../figs/orthocentre.pdf')
plt.savefig('../figs/orthocentre.eps')
#subprocess.run(shlex.split("termux-open ../figs/orthocentre.pdf"))
plt.show()
 

