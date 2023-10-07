% function C = tprod(A,B,transform)
% 
% % C = tprod(A,B,transform) computes the tensor-tensor product of two 3 way tensors A and B under linear transform
% %
% % Input:
% %       A       -   n1*n2*n3 tensor
% %       B       -   n2*m2*n3 tensor
% %   transform   -   a structure which defines the linear transform
% %       transform.L: the linear transform of two types:
% %                  - type I: function handle, i.e., @fft, @dct
% %                  - type II: invertible matrix of size n3*n3
% %
% %       transform.inverseL: the inverse linear transform of transform.L
% %                         - type I: function handle, i.e., @ifft, @idct
% %                         - type II: inverse matrix of transform.L
% %
% %       transform.l: a constant which indicates whether the following property holds for the linear transform or not:
% %                    L'*L=L*L'=l*I, for some l>0.                           
% %                  - transform.l > 0: indicates that the above property holds. Then we set transform.l = l.
% %                  - transform.l < 0: indicates that the above property does not hold. Then we can set transform.l = c, for any constant c < 0.
% %       If not specified, fft is the default transform, i.e.,
% %       transform.L = @fft, transform.l = n3, transform.inverseL = @ifft. 
% %
% %
% % Output: C     -   n1*m2*n3 tensor, C = A * B 
% %
% %
% %
% % See also lineartransform, inverselineartransform
% 
% 
% [n1,n2,n3] = size(A);
% [m1,m2,m3] = size(B);
% 
% if n2 ~= m1 || n3 ~= m3
%     error('Inner tensor dimensions must agree.');
% end
% if nargin < 3
%     % fft is the default transform
%     transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% end
% 
% C = zeros(n1,m2,n3);
% if isequal(transform.L,@fft)
%     % efficient computing for fft transform
% %     A = fft(A,[],3);
% %     B = fft(B,[],3);
%     halfn3 = ceil((n3+1)/2);
%     for i = 1 : halfn3
%         C(:,:,i) = A(:,:,i)*B(:,:,i);
%     end
%     for i = halfn3+1 : n3
%         C(:,:,i) = conj(C(:,:,n3+2-i));
%     end
%     C = ifft(C,[],3);
% else
%     % other transform
%     A = lineartransform(A,transform);
%     B = lineartransform(B,transform);
%     for i = 1 : n3
%         C(:,:,i) = A(:,:,i)*B(:,:,i);
%     end
%     C = inverselineartransform(C,transform);
% end


% function C = tprod(A,B)
% 
% % Tensor-Tensor product of two 3-way tensor: C = A*B
% % A - n1*n2*n3 tensor
% % B - n2*l*n3  tensor
% % C - n1*l*n3  tensor
% %
% % version 1.0 - 18/06/2016
% %
% % Written by Canyi Lu (canyilu@gmail.com)
% % 
% 
% [n1,~,n3] = size(A);
% l = size(B,2);
% A = fft(A,[],3);
% B = fft(B,[],3);
% C = zeros(n1,l,n3);
% for i = 1 : n3
%     C(:,:,i) = A(:,:,i)*B(:,:,i);
% end
% C = real(ifft(C,[],3));
% 


function C=tprod(A,B)
% Computes the tensor-product of two P dimensional tensors (A & B)
% All dimensions of the inputed tensor except for the 2nd must be the same.
% The algorithm follows from paper by ...
%
% INPUTS:
%
% A - n1 x n2 x n3 .... nP tensor
% B - n2 x n4 x n3 .... nP tensor
%
% OUTPUTS: 
%
% C - n1 x n4 x n3 .... nP tensor
%
% Original author :  Misha Kilmer, Ning Hao
% Edited by       :  G. Ely, S. Aeron, Z. Zhang, ECE, Tufts Univ. 03/16/2015


sa = size(A);
sb = size(B);

faces = length(sa);

nfaces = prod(sb(3:end));

for i = 3:faces
    A = fft(A,[],i);
    B = fft(B,[],i);
end
sc = sa;
sc(2) = sb(2);
C = zeros(sc);
for i = 1:nfaces
    C(:,:,i) = A(:,:,i)*B(:,:,i);
end
for i = 3:faces
    C = ifft(C,[],i);
end

