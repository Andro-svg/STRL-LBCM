%==========================================================================
% Compute the tensor contraction of two tensors (Definition 3 in the paper)
% Input:
%   X - the first tensor
%   Y - the second tensor
%   n - The modes of X used to conduct contraction 
%   m - The modes of Y used to conduct contraction 
%  Output:
%   Out - the contraction result
% Created by Yu-Bang Zheng £¨zhengyubang@163.com£©
% Jun. 06, 2020


function Out = tensor_contraction(X,Y,n,m)
Lx = size(X);      Ly = size(Y);   
Nx = ndims(X);     Ny = ndims(Y);
indexx = 1:Nx;     indexy = 1:Ny;
indexx(n) = [];    indexy(m) = [];

tempX = permute(X,[indexx,n]);  tempXX=reshape(tempX,prod(Lx(indexx)),prod(Lx(n)));
tempY = permute(Y,[m,indexy]);  tempYY=reshape(tempY,prod(Ly(m)),prod(Ly(indexy)));
tempOut = tempXX*tempYY;
Out     = reshape(tempOut,[Lx(indexx),Ly(indexy)]);
