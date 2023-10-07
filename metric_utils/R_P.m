function [AUC1,AUC2]=R_P(R,P,w)
%     figure()


x = length(R);
Tao0 = R;
for i = 1:(x-1)
    R(2*i-1)=Tao0(i);
    R(2*i)=Tao0(i);
end
R(2*x-1)=Tao0(x);

y = length(P);
PD_hydice0 = P;
for i = 2:y
    P(2*i-1)=PD_hydice0(i);
    P(2*i-2)=PD_hydice0(i);
end
P(1) = PD_hydice0(1);

AUC1=trapz(R,P);




AUC2 = trapz(R, P);