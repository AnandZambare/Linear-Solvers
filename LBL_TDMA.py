# Line by Line TDMA for any 2-D problem.


import matplotlib.pyplot as plt
from numpy import *

# define TDMA solver here
def TDMA(sub_dia, m_dia, su_dia, rhs):
    # phi is the solution vector.
    v = len(m_dia)
    phi = zeros(v)
    a_s = zeros(v-1)
    b_s = zeros(v-1)
    a_s[0] = m_dia[1] - (sub_dia[1]*su_dia[0]/m_dia[0])
    b_s[0] = rhs[1] - (sub_dia[1]*rhs[0]/m_dia[0])
    # Loop for remaining values.
    for i in range(1, v-1, 1):
        a_s[i] = m_dia[i+1] - (sub_dia[i+1]*su_dia[i]/a_s[i-1])
        b_s[i] = rhs[i+1] - (sub_dia[i+1]*b_s[i-1]/a_s[i-1])
    phi[v-1] = b_s[v-2]/a_s[v-2]
    phi[v-2] = (rhs[v-1]-(m_dia[v-1]*phi[v-1]))/sub_dia[v-1]
    for k in range(v-3, 0, -1):
        phi[k] = (rhs[k+1]-(m_dia[k+1]*phi[k+1])-(su_dia[k+1]*phi[k+2]))/sub_dia[k+1]
    phi[0] = (rhs[0] - (su_dia[0] * phi[1])) / m_dia[0]
    return phi


# define line by line TDMA here which uses TDMA. Still not completely generalized yet. 
def LBLTDMA(n, ny, nx, Mat, R, phi_g):
    phi_n = zeros(n)              # Solution vector
    sub_diagonal = zeros(ny)      # Input for the normal TDMA 
    sub_diagonal[0] = 0.0         # Input for the normal TDMA
    main_diagonal = zeros(ny)     # Input for the normal TDMA
    super_diagonal = zeros(ny)    # Input for the normal TDMA
    super_diagonal[ny-1] = 0.0    # Input for the normal TDMA
    r1 = zeros(ny)                # first factor in the RHS vector.
    r2 = zeros(ny)                # second factor in the RHS vector.
    r3 = zeros(ny)                # third factor in the RHS vector.
    sol = zeros(nx)               # Pseudo vector used to assign the answers.
    # The sweeping directions start from here.
    # Sweeping in Y direction.
    L = 0
    for row in range(0, ny):
        for j in range(0, nx):
            main_diagonal[j] = Mat[L+row+j, L+row+j]
            r1[j] = R[L+row+j]
        for j in range(1, nx):
            sub_diagonal[j] = Mat[L+row+j, row+L+j-1]
        for j in range(0, nx-1):
            super_diagonal[j] = Mat[row+L+j, L+row+j+1]
        # r2 is using the guess values
        if row == ny-1:
        	r2 = zeros(ny)
        else:
            for j in range(0, nx):
            	r2[j] = Mat[L+row+j, L+row+j+nx]*phi_g[L+row+j+nx]  
      
        # r3 is using the newly computed values    
        if row == 0:
        	r3 = zeros(ny)
        else:
            for j in range(0, nx):
            	r3[j] = Mat[L+row+j, L+row+j-nx]*phi_n[L+row+j-nx]
         
        r = r1 - r2 - r3
        sol = TDMA(sub_diagonal, main_diagonal, super_diagonal, r)
        #print(main_diagonal)
        #print(sub_diagonal)
        #print(super_diagonal)
        #print(r)
        #print(sol)
        ## Storing in the new vector.
        for j in range(0, nx):
        	phi_n[j+row+L] = sol[j]
        L = L+(nx-1) ### This the thing to be checked for any random problem. 
        
    # Changing the guess values to the calculated values before changing the direction of sweep.
    for i in range(n):
    	phi_g[i] = phi_n[i]
    # Sweeping in the X direction.
    c1 = zeros(ny)
    c2 = zeros(ny)
    c3 = zeros(ny)
    for col in range(0, ny):
       for j in range(0, nx):
        main_diagonal[j] = Mat[col+nx*j, col+nx*j]
        c1[j] = R[col+nx*j]
       for j in range(1, nx):
        sub_diagonal[j] = Mat[col+nx*j, col+nx*j-nx]
       for j in range(0, nx-1):
        super_diagonal[j] = Mat[nx*j, nx*j+nx+col]
       if col == n-ny:
         c2 = zeros(nx)
       else:
         for j in range(0, nx):
          c2[j] = Mat[nx*j+col, nx*j+col+1]*phi_g[nx*j+col+1]
       if col == 0:
               c3 = zeros(nx)
       else:
         for j in range(0, nx):
              c3[j] = Mat[nx*j+col, nx*j+col-1]*phi_n[nx*j+col-1]
       c = c1-c2-c3
       sol = TDMA(sub_diagonal, main_diagonal, super_diagonal, c)
       for j in range(0, nx):
         phi_n[nx*j+col] = sol[j]
       return phi_n


def Error(sol1, sol2, n):
    res = 0.0
    for i in range(0, n):
            res = res + abs(sol2[i] - sol1[i])
    return res


# Main program.
# Nx = no.of cells in X direction.
# Ny = no.of cells in Y direction.
# N =  no.of cells total in the domain.
# A = Coefficeint Matrix.
# T = transported variable.
# b = RHS or known vector.
N = 16
Nx = 4
Ny = 4
T_g = zeros(N)
A = zeros((N, N))
b = zeros(N)
tol = 1e-4
# Define the coefficients of A here


for i in range  (4, N, 1):
	A[i, i-4] = -1
for i in range(0, N-4, 1):
	A[i, i+4] = -1
for i in range(1, Nx, 1):
	A[i, i-1] = -1
for i in range(13, N, 1):
	A[i, i-1] = -1
for i in range (0, Nx-1, 1):
	A[i, i+1] = -1
for i in range (12, i < N-1, 1):
	A[i, i-1] = -1
for i in range (0, Nx, 1):
	A[i, i] = 5;
A[0, 0] = A[Nx-1, Nx-1] = 6 ;
for i in range(12, N, 1):
	A[i, i] = 5;
A[12, 12] = A[N-1, N-1] = 6 ;
for i in range(Nx, 2*Nx, 1):
	A[i, i] = 4;
A[Nx, Nx] = A[2*Nx-1, 2*Nx-1] = 5 ;
for i in range(2*Nx, 3*Nx, 1):
	A[i, i] = 4;
A[2*Nx, 2*Nx] = A[3*Nx-1, 3*Nx-1] = 5 ;
A[12, 13] = -1;
A[13, 14] = -1 ;
for i in range(Nx, 2*Nx, 1):
	A[i+1][i] = -1
for i in range(2*Nx, 3*Nx, 1):
	A[i+1][i] = -1
for i in range(Nx, 2*Nx-1, 1):
	A[i][i+1] = -1
for i in range(2*Nx, 3*Nx-1, 1):
	A[i][i+1] = -1
A[8, 7] = 0.0 ;
A[12, 11] = 0.0 ;


# Define the b vector here
b[0] =  220 
b[1] = 20  
b[2] = 20 
b[3] = 60  
b[4] = 200 
b[5] = 0 
b[6] = 0 
b[7] = 40 
b[8] = 200 
b[9] = 0 
b[10] = 0 
b[11] = 40 
b[12] =  300 
b[13] = 100 
b[14] = 100  
b[15] = 140  ## RHS of discrete equations
T_n = LBLTDMA(N, Ny, Nx, A, b, T_g)
iteration = 1
Err = Error(T_g, T_n, N)
# print(Err)
while Err > tol:
	# Changing my guess values.
	for i in range(N):
		T_g[i] = T_n[i]
	
	# Solve with new guess values.
	T_n = LBLTDMA(N, Ny, Nx, A, b, T_g)
	# Now calculate the error and continue to work out till solution converges.
	Err = Error(T_g, T_n, N)
	iteration = 1 + iteration
print(T_n)	
print(iteration) 
