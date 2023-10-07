function Out = tnreshape(Grest,N,i)
Nway = size(Grest);
m = zeros(N-1,1);   n = zeros(N-1,1);
for k=1:N-1
    if k<i
        m(k)=2*k;n(k)=2*k-1;
    else
        m(k)=2*k-1;n(k)=2*k;
    end
end
tempG = permute(Grest,[m,n]);
Out = reshape(tempG,prod(Nway(m)),prod(Nway(n)));
