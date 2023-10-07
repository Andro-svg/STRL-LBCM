function Out = tnreshape_new(Grest,N)
Nway = size(Grest);
m = zeros(N-1,1);   n = zeros(N-1,1);
for k=1:N-1
    m(k)=2*k-1;n(k)=2*k;
end
tempG = permute(Grest,[m,n]);
Out = reshape(tempG,prod(Nway(m)),prod(Nway(n)));



