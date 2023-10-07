function [All_Num,time_per_image] = target_detection(readPath, savePath, tuneopts)

if isfield(tuneopts, 'temporal_step');       temporal_step  = tuneopts.temporal_step;   end
if isfield(tuneopts, 'lambdaL');             lambdaL        = tuneopts.lambdaL;         end
if isfield(tuneopts, 'mu');                   mu              = tuneopts.mu;               end


addpath('GoDec_plus\');
filesdir = dir([char(readPath) '/*.jpg']);
if isempty( filesdir )
    filesdir = dir( [char(readPath) '/*.bmp'] );
end
if isempty( filesdir )
    filesdir = dir([char(readPath) '/*.png']);
end

files = { filesdir.name };
files = sort_nat(files);
All_Num = length(files);

iteration = 0;
time_all=0;
S = 0.011;
t_list = [1 : temporal_step : length(files)-temporal_step + 1, length(files)-temporal_step + 1];

t1 = clock;
for t = t_list
    iteration = iteration + 1;
    disp('========================================');
    fprintf('%s %d%s : \n','Starting', iteration, '-th loop');
    spat_temp_ten = [];
    BKG_dict=[];
    BKG_dict_ten = [];
    priorWeight_ten = [];
    BKG = [];   
    tar = [];
    BKG = [];
    Noise = [];
    AO = [];
    r=[];
    %% Read images, obtain prior info and construct image tensor
    for tt=1:temporal_step
        disp([char(readPath) '/' files{tt+t-1}]);
        img = imread([char(readPath) '/' files{tt+t-1}]);
        if size(img, 3) > 1
            img = rgb2gray( img );
        end
        img = double(img);
        imwrite(mat2gray(img), [savePath '/' files{tt+t-1}]);
        [imgHei, imgWid] = size(img);


        %% construct image tensor
        spat_temp_ten(:,:,tt) = img;

        %% Initialize AO for BKG dictionary 
        AO = [AO,img(:)];

    end

    %% Obtain BKG dictionary A tensor
    A = AO';

    [~,Sigma,~]=svd(A,'econ');
    maxsingleV=max(diag(Sigma));
    eta = 0.1;
    r = length(find(diag(Sigma)>maxsingleV*eta));
    q = 0;
    sigma = 1e+4;
    epsilon = 1e-7;
    [L,~,~,~,~] = lowrank_corr(A,r,sigma,epsilon,q,temporal_step);

    T = A - L;
    Num = imgHei*imgWid;
    card=S*Num;

    G = T';
    L = L';

    for tt = 1:size(AO,2)
        BKG_dict = reshape(L(:,tt),[imgHei, imgWid]);
%         imwrite(mat2gray(BKG_dict), [savePath  '/'  strtok([  files{tt+t-1}],'.')  '_dict.jpg']);
        BKG_dict_ten(:, :, tt) = 255*mat2gray(BKG_dict);

        tmp = G(:,tt);
        [~,idx]=sort(abs(tmp(:)),'descend'); 
        tmp(idx(card+1:Num))=0;
        foreground = mat2gray(reshape(tmp, [imgHei, imgWid])); %'GoDec+:foreground')
%         imwrite(foreground, [savePath  '/'  strtok([  files{tt+t-1}],'.')  '_foreground0.jpg']);
        [U,S,V]=svd(foreground,'econ');
        maxS=max(diag(S));
        rr = length(find(diag(S)>maxS*0.001)); %0.005

        [m,~]=size(S);
        S(rr:m,:)=0;
        ref = U*S*V';

        foreground = mat2gray(GuidedFilter(foreground, ref, 1, 0.01));        
%         imwrite(mat2gray(foreground), [savePath  '/'  strtok([  files{tt+t-1}],'.')  '_foregroundguided.jpg']);
        n1=1;
        n2=1;
        k=3*n1+n2+1;
        z = 1.1;
        padded_img = padarray(foreground,[k k],'symmetric');
        prior0 = zeros(imgHei, imgWid);
        for row  = (k+1) : (imgHei + k)
            for col = (k+1) : (imgWid + k) 
                num = 0;
                LBP = '00000000';  

                m0 = padded_img( (row-n1) : (row + n1), (col-n1) : (col + n1));
                mean0 = mean(mean(m0));
         
                m1 = padded_img( (row-3*n1-n2-1) : (row -n1-n2), (col-3*n1-n2-1) : (col -n1-n2));
                m2 = padded_img( (row-3*n1-n2-1) : (row -n1-n2), (col-n1) : (col + n1));
                m3 = padded_img( (row-3*n1-n2-1) : (row -n1-n2), (col+n1+n2) : (col +3*n1+n2+1));
                m4 = padded_img( (row-n1) : (row + n1), (col+n1+n2) : (col +3*n1+n2+1));
                m5 = padded_img( (row+n1+n2) : (row +3*n1+n2+1), (col+n1+n2) : (col +3*n1+n2+1));
                m6 = padded_img( (row+n1+n2) : (row +3*n1+n2+1), (col-n1) : (col +n1));
                m7 = padded_img( (row+n1+n2) : (row +3*n1+n2+1), (col-3*n1-n2-1) : (col -n1-n2));
                m8 = padded_img( (row-n1) : (row + n1), (col-3*n1-n2-1) : (col -n1-n2));

                num = double(mean0>z*max(m1(:)))+double(mean0>z*max(m2(:)))+double(mean0>z*max(m3(:)))+double(mean0>z*max(m4(:)))+double(mean0>z*max(m5(:)))+double(mean0>z*max(m6(:)))+double(mean0>z*max(m7(:)))+double(mean0>z*max(m8(:)));
                
                LBP(8-num+1:8)='1';
                prior0(row-k,col-k) =bin2dec(LBP);     
            end
        end
%         imwrite(mat2gray(prior0), [savePath  '/'  strtok([  files{tt+t-1}],'.')  '_LBCM.jpg']);
        [lambda1, lambda2] = structure_tensor_lambda(prior0, 'Gaussian', 3);
        x=3;
        y=1;
        prior =(x+y)*lambda1.*lambda2./(y*lambda2+x*lambda1+0.001); % corner extract
%         imwrite(mat2gray(prior), [savePath   '/'  strtok([  files{tt+t-1}],'.')  '_prior.jpg']);
        priorWeight_ten(:,:,tt) = prior;
    end
    
    %% The LRSD model
    Nway = size(spat_temp_ten);
    N = ndims(spat_temp_ten);
    
    % Paramter Setting
    opts=[];
    opts.max_iter = 400;
    opts.tol =1e-4;
    omega = cell(1,N);
    weight = 4;
    omega{1} = 1/(2+weight);
    omega{2} = 1/(2+weight);
    omega{3} = weight/(2+weight);
    opts.omega = omega;
    opts.max_R = 20 * [0,  1,  1;
                       0,  0,  1;
                       0,  0,  0];
    opts.R     = 3 * [0,  1,  1;
                      0,  0,  1;
                      0,  0,  0];

    opts.gamma = 1.4;
    opts.lambda = 1;
    opts.lambda1 = lambdaL*1000/(sqrt(min([Nway(1),Nway(2),Nway(3)])));
    opts.lambda2 = 50* opts.lambda1;
    opts.lambda3 =  0.1*opts.lambda1;
    opts.beta = 1;
    opts.alpha = 0.5;
    opts.mu  = mu;
    opts.max_mu = 1e10;
    tenT=[];
    tenB=[];
    tenN=[];
    [tenB, tenT, tenN] = LRSD(spat_temp_ten, BKG_dict_ten, priorWeight_ten, opts);

    %% reconstrcut imgNum images
    for kk = 1:temporal_step  
        tar = tenT(:,:,kk);
        BKG = tenB(:,:,kk);
        Noise = tenN(:,:,kk);
        E = tar;
        imwrite(mat2gray(E), [savePath  '/' strtok([files{kk+t-1}],'.') '_tar.jpg']);
%         imwrite(mat2gray(BKG), [savePath '/result_' num2str(C) '/' strtok([files{kk+t-1}],'.') '_BKG.jpg']);
%         imwrite(mat2gray(Noise), [savePath '/result_' num2str(C) '/' strtok([files{kk+t-1}],'.') '_Noise.jpg']);
    end

end
t2=clock;
disp(['Programming has running:',num2str(etime(t2,t1))]);
disp('=====================================================')
time_all = etime(t2,t1);
time_per_iter = time_all/iteration;
disp(['Each iteration consumes time: ', num2str(time_per_iter)]);
time_per_image = time_per_iter / temporal_step;
disp(['Each image consumes time: ', num2str(time_per_image)]);
end