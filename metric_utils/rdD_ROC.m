function [t]=rdD_ROC(Tao,PF_hydice,PD_hydice,w)
%     figure()
%     plot3(PF_hydice(:),Tao(:),PD_hydice(:),'-r','LineWidth',1);
%     axis([0 1 0 1 0 1])
%     grid on
%     title ('3D ROC')
%     xlabel('P_{F}');ylabel('\tau');zlabel('P_{D}');
%   legend ('RX-AD','R-AD');
t=1;