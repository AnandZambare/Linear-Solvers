#include <iostream>
#include <math.h>
#include <vector>
using namespace std;
void printMatrix (double** mat, int m) {
for (int i = 0; i < m; i ++) {
for (int j = 0; j < m; j++) {
cout<< mat[i][j]<<" ";
}
cout <<endl;
}
return;
}
void printVector (vector<double> vec){
for (int j=0; j < vec.size(); j++){
cout<< vec[j]<<" " ;
}
return ;
}
double error(vector<double> Z_old, vector<double> Z_new){
double res = 0.0 ;
for (int i = 0; i < Z_old.size() ; i++){
res = res + abs( Z_new[i] - Z_old[i] ) ;
}
return res ;
}
vector<double> old_fac(double** mat, vector<double> x1){
// this calculates the old factor
int q = x1.size() ;
int m, n ;
vector<double> fac_1(q) ;
for (m = 0; m < q; m++){
fac_1[m] = 0.0 ;
for (n = m+1; n < q; n++){
fac_1[m] = (fac_1[m] + (mat[m][n]*x1[n])) ;
}fac_1[m] = fac_1[m]/mat[m][m] ;
}
fac_1[q-1] = 0.0 ;
return fac_1 ;
}
int main(){
int N, i, j, k, p;
cout << "Enter the no.of equations"<<endl;
cin>>N;
vector<double> b(N), phi_g(N), phi_n(N), r1(N), r2(N), r3(N);
// defining the coefficient matrix
double** A = new double*[N];
for (i = 0; i < N; i++) {
A[i] = new double[N];
}
// specific for 3 by 3 system
b[0] =  351.84126984; b[1] = 120.15873016 ; b[2] = 154.50793651;    // RHS of discrete equations
// Coefficeint Matrix
A[0][0] = 6;  A[0][1] = -1; A[0][2] = 0;  
A[1][0] = -1; A[1][1] = 5;  A[1][2] = -1; 
A[2][0] = 0;  A[2][1] = -1; A[2][2] = 6;  
cout << endl ;
// printMatrix(A, N);
// set the tolerance value and also set the guess value.
double tol = 1e-3, err;
for (i = 0; i < N; i++){
phi_g[i] = 0.0 ;
}
// calculate the constant part of the rhs in solved equations.
for (i = 0; i < N; i++){
r1[i] = b[i]/A[i][i] ;
}// calculate the old factor part.
r2 = old_fac(A, phi_g);
// r3 is the dynamic part should be evaluated with the new vector only.
r3[0] = 0.0 ;
phi_n[0] = r1[0] - r2[0] - r3[0] ;
// dynamic calculations for r3 and phi_n
for (k = 1; k < N; k++){
r3[k] = 0.0 ;
for (j = 0; j < k+1 ; j++){
r3[k] = (r3[k] + (A[k][j]*phi_n[j])) ;
}
r3[k] = r3[k]/A[k][k] ;
phi_n[k] = r1[k] - r2[k] - r3[k] ;
}
err = error(phi_g, phi_n);
// now go throught the loop.
while (err > tol){
for (i = 0; i < N; i++){
phi_g[i] = phi_n[i] ;
phi_n[i] = 0.0 ;
}
// now update the old factors
r2 = old_fac(A, phi_g) ;
r3[0] = 0.0 ;
phi_n[0] = r1[0] - r2[0] - r3[0] ;
for (k = 1; k < N; k++){
r3[k] = 0.0;
for (j = 0; j < k+1; j++){
r3[k] = (r3[k] + (A[k][j]*phi_n[j])) ;
}
r3[k] = r3[k]/A[k][k] ;
phi_n[k] = r1[k] - r2[k] - r3[k] ;
}
err = error(phi_g, phi_n);
}
//printVector(r3);
//cout<<endl ;
printVector(phi_n);
cout<<endl;
cout<< err ;
cout<< endl ;
return 0 ;
}
