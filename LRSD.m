function [B,T,Noise] = LRSD(D, A, P, opts)

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');            tol      = opts.tol;                  end
if isfield(opts, 'max_iter');       max_iter = opts.max_iter;             end
if isfield(opts, 'lambda');         lambda   = opts.lambda;               end
if isfield(opts, 'lambda1');        lambda1  = opts.lambda1;              end
if isfield(opts, 'lambda2');        lambda2  = opts.lambda2;              end
if isfield(opts, 'lambda3');        lambda3  = opts.lambda3;              end
if isfield(opts, 'omega');          omega    = opts.omega;                end
if isfield(opts, 'alpha');          alpha    = opts.alpha;                end
if isfield(opts, 'beta');           beta     = opts.beta;                 end
if isfield(opts, 'R');              R        = opts.R;                    end
if isfield(opts, 'max_R');          max_R    = opts.max_R;                end
if isfield(opts, 'gamma');          gamma    = opts.gamma;                end
if isfield(opts, 'mu');             mu       = opts.mu;                   end
if isfield(opts, 'max_mu');         max_mu   = opts.max_mu;               end
 
%% Initialization
Nway = size(D);
N = ndims(D);
Noise = zeros(Nway);
B = ones(Nway);
Xway = [Nway(2),Nway(2),Nway(3)];
tempdim = diag(Xway)+R+R';
max_tempdim = diag(Nway)+max_R+max_R';

G = cell(1,N);
Q = cell(1,N);
for i = 1:N
    G{i} = ones(tempdim(i,:));
    Q{i} = ones(tempdim(i,:));
end

T = zeros(Nway);
V1 = zeros(Xway);
V2 = ones(Nway);
V3 = ones(Nway);
V4 = ones(Nway);
M1 = zeros(Nway);
M2 = zeros(Xway);
M3 = zeros(Nway);
M4 = zeros(Nway);
M5 = zeros(Nway);
Z = cell(1,N);
for i=1:N
    Z{i} = zeros(size(G{i})); % M_{6,k}(k=1,2,...,N)
end
M7 = zeros(Nway);

preNumT = numel(T);
W = ones(size(T))./ (P+0.001);

m1 = [1 0;0 0 ];
m2 = [-1 0 ;0 0 ];
template_time(:,:,1) = m1;
template_time(:,:,2) = m2;

% horizontal difference operators
FDx = psf2otf([1 -1],Nway);
% vertical difference operator
FDy = psf2otf([1;-1],Nway);
% z difference operator
FDz = psf2otf(template_time,Nway);

FDxH = conj(FDx);
FDyH = conj(FDy);
FDzH = conj(FDz);

AT = tprod(tran(A),A);
I = teye(Nway(2),Nway(3));
IL = 1./((abs(FDx).^2 + abs(FDy).^2 + abs(FDz).^2) + 1);

for iter = 1 : max_iter
    disp(' ')
    %% update X
    temp1 = tprod(tran(A),(D-T-Noise-M1+B-M7))+(V1-M2);
    inverse = t_inverse(I + 2 * AT);
    X = tprod(inverse,temp1);
    AX = tprod(A,X);
    
    %% update Gk
    for i = 1:N
        V1i = my_Unfold(V1,size(V1),i);
        Qii = my_Unfold(Q{i},size(Q{i}),i);
        Zii = my_Unfold(Z{i},size(Z{i}),i);
        Gi = my_Unfold(G{i},tempdim(i,:),i);
        Girest = tnreshape(tnprod_rest(G,i),N,i);
        tempC = mu*(Qii-Zii)+lambda*V1i*Girest';
        tempA = lambda*(Girest*Girest')+mu*eye(size(Gi,2));
        G{i}  = real(my_Fold(tempC*pinv(tempA),tempdim(i,:),i));
    end

    %% update Qk
    for i = 1:N
        [Q{i},~,~] = prox_tnn_my(G{i}+Z{i},omega{i}/mu);
        Q{i} = real(Q{i});
    end

    %% update B
    temp_Bk = AX+M7;
    FB = IL.*(FDxH.*fftn(V2-M3) + FDyH.*fftn(V3-M4) + FDzH.*fftn(V4-M5) + fftn(temp_Bk));
    B = real(ifftn(FB));

    %% update V1
    V1 = (lambda*tnprod(G)+mu*(X+M2))/(lambda+mu); 

    %% update V2 and V3
    tempV2 = ifftn(FDx.*FB)+M3+lambda3*beta*V2./(mu*sqrt(V2(:)' * V2(:) + V3(:)' * V3(:)));
    V2 = prox_l1(tempV2,lambda3/mu);
  
    tempV3 = ifftn(FDy.*FB)+M4+lambda3*beta*V3./(mu*sqrt(V2(:)' * V2(:) + V3(:)' * V3(:)));
    V3 = prox_l1(tempV3,lambda3/mu);
    
    %% update V4
    tempV4 = ifftn(FDz.*FB)+M5;
    V4 = prox_l1(tempV4,alpha*lambda3/mu);

    %% update T
    tempT = D-AX-Noise-M1;
    thres = W.*lambda1/mu;
    T = real(prox_l1(tempT, thres));
    
    %% update W
    W = 2 ./ ((abs(T))+ 0.01) ./ (P+0.001);

    %% update Noise 
    Noise = real((mu*(D-AX-T-M1))/(2*lambda2+mu));

    %% check the convergence
    currNumT = sum(T(:) > 0); 
    chg =norm(D(:)-B(:)-T(:)-Noise(:))/norm(D(:));
    fprintf('iter = %d   res=%.10f  \n', iter, chg);    

    if (chg < tol) || (currNumT == preNumT)
        break;
    end
    preNumT = currNumT;
    
    %% update Lagrange multipliers M and penalty parameter mu
    M1 = M1 - (D-AX-T-Noise);
    M2 = M2 - (V1-X);
    M3 = M3 - (V2-ifftn(FDx.*FB));
    M4 = M4 - (V3-ifftn(FDy.*FB));
    M5 = M5 - (V4-ifftn(FDz.*FB));
    M7 = M7 - (B-AX);
    for i=1:N
        Z{i} = Z{i} - (Q{i}-G{i});
    end
    mu = min(gamma*mu,max_mu);  

    %% update the estimated rank
    rank_inc=double(tempdim<max_tempdim);
    if sum(rank_inc(:))~=0
        G = rank_inc_adaptive(G,rank_inc,N);
        Q = rank_inc_adaptive(Q,rank_inc,N);
        Z = rank_inc_adaptive(Z,rank_inc,N);
        tempdim = tempdim+rank_inc;
    end
end

end

function [G]=rank_inc_adaptive(G,rank_inc,N)
    % increase the estimated rank
    for j = 1:N
        G{j} = padarray(G{j},rank_inc(j,:),1,'post');
    end
end
