# This is the Tridiagonal matrix solver.

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
    
    
# Just provide the three diagonals to the solver with Right hand side vector in the main code and get the direct solution.
