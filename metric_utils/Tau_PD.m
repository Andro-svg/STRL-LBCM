function [AUC]=Tau_PD(Tao,PD_hydice,w)
%     figure()


% x = length(Tao);
% Tao0 = Tao;
% for i = 2:x
%     Tao(2*i-1)=Tao0(i);
%     Tao(2*i-2)=Tao0(i);
% end
x = length(Tao);
Tao0 = Tao;
for i = 1:(x-1)
    Tao(2*i-1)=Tao0(i);
    Tao(2*i)=Tao0(i);
end
Tao(2*x-1)=Tao0(x);


% y = length(PD_hydice);
% PD_hydice0 = PD_hydice;
% for i = 1:(y-1)
%     PD_hydice(2*i-1)=PD_hydice0(i);
%     PD_hydice(2*i)=PD_hydice0(i);
% end
% PD_hydice(2*y-1) = PD_hydice0(y);
% AUC=trapz(Tao,PD_hydice);

y = length(PD_hydice);
PD_hydice0 = PD_hydice;
for i = 2:y
    PD_hydice(2*i-1)=PD_hydice0(i);
    PD_hydice(2*i-2)=PD_hydice0(i);
end
PD_hydice(1) = PD_hydice0(1);

% *
% |
% |_ _ _ *
%        |
%        |
%        |_ _ _ *
AUC=trapz(Tao,PD_hydice);

%     plot(Tao(:),PD_hydice(:),'-r','LineWidth',1);
%     axis([0 1 0 1])
%     grid on
%     title ('P_{D} vs \tau')
%     xlabel('\tau');ylabel('P_{D}');
%   legend ('RX-AD','R-AD');

% for j=1:w 
%     AUC0=trapz(Tao(:,j),PD_hydice(:,j));
%     a=PD_hydice(size(PD_hydice,1),j);
%     AUC(j,1)=(AUC0-a)/(1-a);
% end