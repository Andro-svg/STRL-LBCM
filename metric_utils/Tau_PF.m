function [AUC]=Tau_PF(Tao,PF_hydice,w)
%     figure()


x = length(Tao);
Tao0 = Tao;
for i = 1:(x-1)
    Tao(2*i-1)=Tao0(i);
    Tao(2*i)=Tao0(i);
end
Tao(2*x-1)=Tao0(x);

y = length(PF_hydice);
PF_hydice0 = PF_hydice;
% for i = 1:(y-1)
%     PF_hydice(2*i-1)=PF_hydice0(i);
%     PF_hydice(2*i)=PF_hydice0(i);
% end
% PF_hydice(2*y-1) = PF_hydice0(y);
for i = 2:y
    PF_hydice(2*i-1)=PF_hydice0(i);
    PF_hydice(2*i-2)=PF_hydice0(i);
end
PF_hydice(1) = PF_hydice0(1);

% *
% |
% |_ _ _ *
%        |
%        |
%        |_ _ _ *

AUC=trapz(Tao,PF_hydice);
    
%     plot(Tao(:),PF_hydice(:),'-r','LineWidth',1);
%     %stairs(Tao(:),PF_hydice(:),'-r','LineWidth',1);
%     axis([0 1 0 1])
%     grid on
%     title ('P_{F} vs \tau')
%     xlabel('\tau');ylabel('P_{F}');
%   legend ('RX-AD','R-AD');
% for j=1:w
%     AUC0=trapz(Tao(:,j),PF_hydice(:,j));
%     b=PF_hydice(1,j);
%     AUC(j,1)=AUC0/b;
% end
    