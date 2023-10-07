function [Precision, Recall ,F1,AUC1,AUC2]=pr_curve(GTpath, ImgPath, seq_name,method_name,ImgNum, Slack, thres_mode)

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

%% 获取该数据集检测结果的TP像素个数GT
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
R=zeros(length(Thres),1);
P=zeros(length(Thres),1);
%用于记录总的真实目标个数和总的检测到的目标的个数，总的检测到的目标的个数：所有tau下，所有图像中检测到的目标的总数
true_target_sum = 0; % tau的个数*（所有图像的真实目标个数总和）
detected_target_sum = 0; % 每个tau下面得到的所有图像的真实检测到的目标个数总和
BKG_pixel_sum =0;
false_detect_sum = 0;


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
        offset_pos = offset_pos + length(gt_x_mids); % 每一张图的真实目标个数相加，最终获得真实的总目标个数offset_pos，用于计算TP
        false_detect = false_detect + false_detect_per; % 每一张图的误检像素数增加，最终获得总错误检测像素数量， 用于计算FP
    end
    FPR(i)=false_detect/(divisor-TP_pixel);  % 错误检测的像素点总数/背景像元个数-TP_pixel    divisor=divisor+I_N*I_M是数据集图像总面积   
    TPR(i)=sum(TP)/offset_pos;  % sum(TP)一个tau下面得到的检测正确的目标的个数，    真实的总目标个数offset_pos


    P(i)=sum(TP)/(false_detect+sum(TP));
    if isnan(P(i))
        P(i)=1;
    end
    R(i)=sum(TP)/offset_pos;
   
end

P=[1;P;0];
R=[0;R;1];
Tau_pr=flip([1,Thres,0]');


TP_pr=zeros(10000,1);
offset_pos_pr=0;
false_detect_pr=0;
for img_ind=1:ImgNum
%         gt_arr=read_gt_xml(P,GTpath, seq_name, img_ind);  % gt_arr(i,:)=[xmin,ymin,xmax,ymax];
    gt_arr=read_gt_txt(GTpath, img_ind, img_name);
%         gt_arr=read_gt_txt_v2(GTpath, seq_name, img_ind);
    TheMap=read_detect_map(ImgPath, seq_name,img_ind, method_name);  
    M = mean(TheMap(:));
    sigma = std(TheMap(:));
    thres = M + sigma;

    gt_x_mids=(gt_arr(:,1)+gt_arr(:,3))/2;
    gt_y_mids=(gt_arr(:,2)+gt_arr(:,4))/2;
    [Trow,Tcol]=find(TheMap>=thres);  % y,x
    %%%%%%%%%%%%%%更改于2023-09-14%%%%%%%%%%%%%%%
    % 像素点落在以目标为中心，Slack尺度范围内算正确检测，其余都算误检
    false_detect_per_pr = length(Trow);  % 当前阈值下切割得到的第img_ind张图像的像素数量
    for jj=1:length(gt_x_mids)  % 逐个ground_truth与每个像素点进行判断
        for ii=1:length(Trow)  % 每个像素点与第jj个gt进行距离判断
            Tdis=max(abs(Trow(ii)-gt_y_mids(jj)),abs(Tcol(ii)-gt_x_mids(jj)));
            if Tdis<=Slack
                TP_pr(offset_pos_pr+jj)=1;  %如果有像素点落在第jj个gt的slack范围内，则第jj个目标点检测正确
                false_detect_per_pr = false_detect_per_pr-1; % 错误检测的像素点数量-1
            end                
        end
    end
    offset_pos_pr = offset_pos_pr + length(gt_x_mids); % 每一张图的真实目标个数相加，最终获得真实的总目标个数offset_pos，用于计算TP
    false_detect_pr = false_detect_pr + false_detect_per_pr; % 每一张图的误检像素数增加，最终获得总错误检测像素数量， 用于计算FP
end

% precision
Precision=sum(TP_pr)/(false_detect_pr+sum(TP_pr));
% recall
Recall=sum(TP_pr)/offset_pos_pr;
% F1
F1=2*Precision*Recall./(Recall+Precision);

[AUC1,AUC2]=R_P(R,P,1);

save([ImgPath '/P.mat'],'P');
save([ImgPath '/R.mat'],'R');
save([ImgPath '/Tau_pr.mat'],'Tau_pr');

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


function gt_arr=read_gt_txt_v2(GTpath,Seq_name, img_ind)
readPath = GTpath;
filesdir = dir([GTpath '/*.txt']);
files = { filesdir.name };
files = sort_nat(files);
imagesdir = dir([GTpath '/*.png']);
images = { imagesdir.name };
img = imread([GTpath '/' images{1}]);
[h,w] = size(img);

gt_pos= load([readPath '/' char(strtok(files(img_ind),'.')) '.txt']);
[tar_num,~] = size(gt_pos);
gt_arr=zeros(tar_num,4);
for i=1:tar_num
    xmin= w * gt_pos(i,2);
    ymin= h * gt_pos(i,3);
    xmax= w * gt_pos(i,2);
    ymax= h * gt_pos(i,3);
    gt_arr(i,:)=[xmin,ymin,xmax,ymax];
end
end


function TheMap=read_detect_map(ImgPath, seq_name, img_ind, method_name)
% if strcmp(method_name{1},'IVSTTM')
%     filesdir = dir([ImgPath '/*' char(method_name) '_target.mat']);
% else
%     filesdir = dir([ImgPath '/*' char(method_name) '.mat']);
% end
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
