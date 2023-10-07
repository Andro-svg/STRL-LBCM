function  Xout = end_permute(X,k)
N     = ndims(X);
n     = 1:N;
n(k)  = [];
index = [n,k]; 
Xout  = permute(X,index);