%==========================================================================
% Compute the FCTN composition of the factors (Definition 5 in the paper)
% Input:
%   G - the FCTN factors 
%  Output:
%   Out - the FCTN composition of G
% Created by Yu-Bang Zheng £¨zhengyubang@163.com£©
% Jun. 06, 2020

function Out = tnprod_new(X)
N = length(X);
m = 1; n = 1;
Out = X{1};
for i = 1:N-1
    Out = tensor_contraction(Out,X{i+1},m,n);
    n = [n,1+i];
    tempm = 1+i*(N-i);
    if i>1
        m(2:end) = m(2:end)-[1:i-1];
    end
    m   = [m, tempm];
end