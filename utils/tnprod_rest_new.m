%==========================================================================
% Compute the FCTN composition of N-1 FCTN factors, i.e., one factor 
% does not participate in the composition (Definition 5 in the paper)
% Input:
%   G - the FCTN factors 
%   rest - the index of the factor that does not participate in the composition
%  Output:
%   Out - the FCTN composition of G_1,G_2,...,G_rest-1,G_rest+1,...,G_N
% Created by Yu-Bang Zheng £¨zhengyubang@163.com£©
% Jun. 06, 2020

function Out = tnprod_rest_new(G,rest)
N = length(G);
m1 = zeros(1,N-2); m2 = zeros(1,N-2);
n = 1:N-1; 
if rest < N
   n(rest)=[];
end
for i=1:N-2
    m1(i)=1+(i-1)*N;
    m2(i)=2+(i-1)*N;
end
j=1;
if rest>1
    Out = G{1};
    for i=1:N-1
        if i+1 < rest
            Out = tensor_contraction(Out,G{i+1},m1(1:j),n(1:j));
            m1(2:end)=m1(2:end)-[1:N-3];
            m2(2:end)=m2(2:end)-[1:N-3];
            j=j+1;
        end
        if i+1 > rest
            Out = tensor_contraction(Out,G{i+1},m2(1:j),n(1:j));
            m2(2:end)=m2(2:end)-[1:N-3];
            j=j+1;
        end
    end
end
if rest == 1
    Out = G{2};
    for i=2:N-1
        Out = tensor_contraction(Out,G{i+1},m2(1:i-1),n(1:i-1));
        m2(2:end)=m2(2:end)-[1:N-3];
    end
end