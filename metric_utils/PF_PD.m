function [AUC]=PF_PD(PF_hydice,PD_hydice,w)
%     figure()
    

% x = length(PF_hydice);
% PF_hydice0 = PF_hydice;
% for i = 2:x
%     PF_hydice(2*i-1)=PF_hydice0(i);
%     PF_hydice(2*i-2)=PF_hydice0(i);
% end

x = length(PF_hydice);
PF_hydice0 = PF_hydice;
for i = 1:(x-1)
    PF_hydice(2*i-1)=PF_hydice0(i);
    PF_hydice(2*i)=PF_hydice0(i);
end
PF_hydice(2*x-1)=PF_hydice0(x);


% y = length(PD_hydice);
% PD_hydice0 = PD_hydice;
% for i = 1:(y-1)
%     PD_hydice(2*i-1)=PD_hydice0(i);
%     PD_hydice(2*i)=PD_hydice0(i);
% end
% PD_hydice(2*y-1) = PD_hydice0(y);
y = length(PD_hydice);
PD_hydice0 = PD_hydice;
for i = 2:y
    PD_hydice(2*i-1)=PD_hydice0(i);
    PD_hydice(2*i-2)=PD_hydice0(i);
end
PD_hydice(1) = PD_hydice0(1);
%                 _ _ _ *
%                |
%          _ _ _ * 
%        |
%        |
%        * 
AUC=abs(trapz(PF_hydice,PD_hydice));

%     plot(PF_hydice(:),PD_hydice(:),'-r','LineWidth',1); 
%     axis([0 1 0 1])
%     grid on
%     title ('P_{D} vs P_{F}')
%     xlabel('P_{F}');ylabel('P_{D}');
%   legend ('RX-AD','R-AD');
% for j=1:w
%         AUC0=abs(trapz(PF_hydice(:,j),PD_hydice(:,j)));
%         a=PD_hydice(size(PD_hydice,1),j);
%         AUC(j,1)=(AUC0-a)/(1-a);        
% end