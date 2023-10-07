function [Tau,FPR,TPR,AUC_pfpd,AUC_tpd,AUC_tpf,ODP,TD,BS,TDBS,SNPR]=three_d_roc_curve(GTpath, ImgPath, seq_name,method_name,ImgNum, Slack, thres_mode)
%{
统一用同一个阈值对同一序列所有map进行分割计算PR
thres_mode=1时，阈值采用每张检测结果的像元的pixel value
thres_mode=2时，阈值采用等差数列
parameters：
--seq_name：string format, sequence name
--method_name: string format, method name
--ImgeNum: double or uint8 format, number of sequence images
--Slack: double or uint8 format, a slack size for target region, default as 0  计算检测到的点到中心的切比雪夫距离，如果距离大于Slack，是虚警
--thres_mode: uint8 format, threshold definition mode
%}
if nargin<6
    thres_mode=1;
end

%% 获取阈值
divisor=0;
Thres=[];
if thres_mode==1
    for img_ind=1:ImgNum
        TheMap=read_detect_map(ImgPath, seq_name,img_ind, method_name);
        [I_N,I_M]=size(TheMap);
%         max(max(TheMap))
        divisor=divisor+I_N*I_M;
        tmp = reshape(TheMap,I_N*I_M,1);
%         tmp=sort(tmp,'descend');
        Thres=[Thres;tmp];
    end
    Thres=unique(Thres);
else
    for img_ind=1:ImgNum
        TheMap=read_detect_map(ImgPath, seq_name,img_ind, method_name);
        [I_N,I_M]=size(TheMap);
        divisor=divisor+I_N*I_M;   % 数据集图像总面积
    end
%    Thres=round((1:-0.01:0.01)*255);
    Thres=1:-0.004:0;
end
Thres=sort(Thres,'descend');
                                                                                                                                                                                             
%% 获取图像名称列表
filesdir = dir([char(GTpath) '/*.jpg']);
if isempty( filesdir )
    filesdir = dir( [char(GTpath) '/*.bmp'] );
end
if isempty( filesdir )
    filesdir = dir([char(GTpath) '/*.png']);
end
files = { filesdir.name };
files = sort_nat(files);
img_name=[];
for kk=1:length(files)
    temp_name = files(kk);
    temp_name = temp_name{1};
    img_name=[img_name str2num(strtok(temp_name,'.'))];
end

%% 获取该数据集的TP像素个数
TP_pixel = 0;
for img_ind=1:ImgNum  %一共有ImgNum张图片,第img_ind张图片
%     gt_arr=read_gt_xml(P,GTpath, seq_name, img_ind);  % gt_arr(i,:)=[xmin,ymin,xmax,ymax];
    gt_arr=read_gt_txt(GTpath, img_ind, img_name);
%     gt_arr=read_gt_txt_v2(GTpath, seq_name, img_ind);
    TheMap=read_detect_map(ImgPath, seq_name,img_ind, method_name);
    [I_N,I_M]=size(TheMap);
    template0 = ones(I_N,I_M);
    template = padarray(template0,[Slack Slack],0,'both');
    gt_x_mids=(gt_arr(:,1)+gt_arr(:,3))/2;
    gt_y_mids=(gt_arr(:,2)+gt_arr(:,4))/2;
    for jj=1:length(gt_x_mids)
        tmp_x = gt_x_mids(jj)+Slack;
        tmp_y = gt_y_mids(jj)+Slack;
        tmp = template((tmp_x-Slack):(tmp_x+Slack),(tmp_y-Slack):(tmp_y+Slack));
        TP_pixel = TP_pixel + sum(tmp(:));
    end
end


%% 计算TPR和FPR
TPR=zeros(length(Thres),1);
FPR=zeros(length(Thres),1);
for i=1:length(Thres)
    thres=Thres(i);
    TP=zeros(10000,1);
    offset_pos=0;
    false_detect=0;
    for img_ind=1:ImgNum
%         gt_arr=read_gt_xml(P,GTpath, seq_name, img_ind);  % gt_arr(i,:)=[xmin,ymin,xmax,ymax];
        gt_arr=read_gt_txt(GTpath, img_ind, img_name);
%         gt_arr=read_gt_txt_v2(GTpath, seq_name, img_ind);
        TheMap=read_detect_map(ImgPath, seq_name,img_ind, method_name);                    
        gt_x_mids=(gt_arr(:,1)+gt_arr(:,3))/2;
        gt_y_mids=(gt_arr(:,2)+gt_arr(:,4))/2;
        [Trow,Tcol]=find(TheMap>=thres);  % y,x
        %%%%%%%%%%%%%%更改于2023-01-17%%%%%%%%%%%%%%%
        % 像素点落在以目标为中心，Slack尺度范围内算正确检测，其余都算误检
        false_detect_per = length(Trow);  % 当前阈值下切割得到的第img_ind张图像的像素数量
        for jj=1:length(gt_x_mids)  % 逐个ground_truth与每个像素点进行判断
            for ii=1:length(Trow)  % 每个像素点与第jj个gt进行距离判断
                Tdis=max(abs(Trow(ii)-gt_y_mids(jj)),abs(Tcol(ii)-gt_x_mids(jj)));
                if Tdis<=Slack
                    TP(offset_pos+jj)=1;  %如果有像素点落在第jj个gt的slack范围内，则第jj个目标点检测正确
                    false_detect_per = false_detect_per-1; % 错误检测的像素点数量-1
                end
            end
        end
        offset_pos = offset_pos + length(gt_x_mids); % 每一张图的目标个数相加，最终获得总目标个数，用于计算TP
        false_detect = false_detect + false_detect_per; % 每一张图的误检像素数增加，最终获得总错误检测像素数量， 用于计算FP
    end
    FPR(i)=false_detect/(divisor-TP_pixel);  % 错误检测的像素点总数/背景像元个数-TP_pixel             
    TPR(i)=sum(TP)/offset_pos;
end


%% ***********************************3D ROC
Tau=flip(Thres');
FPR=flip(FPR);TPR=flip(TPR);
[AUC_pfpd]=PF_PD(FPR,TPR,1);
[AUC_tpd]=Tau_PD(Tau,TPR,1);
[AUC_tpf]=Tau_PF(Tau,FPR,1);            
[~]=rdD_ROC(Tau,FPR,TPR,1);
save([ImgPath '/FPR.mat'],'FPR');
save([ImgPath '/TPR.mat'],'TPR');
save([ImgPath '/Tau.mat'],'Tau');

ODP=AUC_pfpd+AUC_tpd-AUC_tpf;
TD=AUC_pfpd+AUC_tpd;
BS=AUC_pfpd-AUC_tpf;
TDBS=AUC_tpd-AUC_tpf;
SNPR=AUC_tpd./AUC_tpf;
end


function gt_arr=read_gt_txt(GTpath, img_ind, img_name)
readPath = GTpath;
filesdir = dir([GTpath '/*.txt']);
files = { filesdir.name };
gt_pos = load([readPath '/' char(files)]);
tar_num = 1;
gt_arr=zeros(tar_num,4);
for i=1:tar_num
    xmin=gt_pos(img_name(img_ind)+1,2);
    ymin=gt_pos(img_name(img_ind)+1,3);
    xmax=gt_pos(img_name(img_ind)+1,2);
    ymax=gt_pos(img_name(img_ind)+1,3);
    gt_arr(i,:)=[xmin,ymin,xmax,ymax];
end
end


function TheMap=read_detect_map(ImgPath, seq_name, img_ind, method_name)
filesdir = dir([ImgPath '/*' char(method_name) '.mat']);
files = { filesdir.name };
files = sort_nat(files);
mapname=fullfile(ImgPath,'/',strcat(strtok(files(img_ind),'.'),'.mat'));
TheMap0=load(char(mapname));
% TheMap = mat2gray(TheMap0.E);
% TheMap=round((TheMap-min(TheMap(:)))/(max(TheMap(:))-min(TheMap(:))+eps)*255);  %归一化为8位图像
TheMap = TheMap0.E;
TheMap = double(TheMap);
TheMap=(TheMap-min(TheMap(:)))/(max(TheMap(:))-min(TheMap(:))+eps); % eps防止分母为零
end
