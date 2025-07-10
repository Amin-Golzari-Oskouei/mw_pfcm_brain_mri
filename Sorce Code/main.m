% This demo shows how to call the Neighborhood based Feature & Cluster Weighted PFCM algorithm described in the paper:
% A.Najafi, N.Aghazadeh, M.Hashemzadeh and A.Golzari oskouei,
% "An Adaptive Multi-Weighted Possibilistic Fuzzy C-Means Clustering Approach for Brain Magnetic Resonance Imaging Segmentation"
function [accurcy,fm,nmi, Jaccard, Sensitivity] = main(FileName, PathName)
%%
f_ori=imread(strcat(PathName,'\original','\',FileName));
GT = imread(strcat(PathName,'\GT','\',FileName));
% figure,imshow(GT);
%% skull border remove
binarymak = skull_remove(f_ori);
f_ori = f_ori(3:end-3, 4:end-4, :);
GT = GT (3:end-3, 4:end-4);
f_ori = bsxfun(@times, f_ori, cast(binarymak, 'like', f_ori));
%% Curvelet Transform
f_ori=EdgeEnhancedCurvelet(f_ori);
f_ori=bilatFilter(f_ori);
%  figure,imshow(f_ori);
%% Feathre Extract step
fprintf('The feature extraction phase has started ...\n')
X = FeatureExtractor(f_ori);
X(:, 1:3) = [];
[N,d]=size(X);
X=(X(:,:)-min(X(:)))./(max(X(:)-min(X(:)))); %Normalize data between 0 and 1 (optinal)

%Algorithm parameters.
%---------------------
k=3;  %number of clusters.
q=-8;                     %the value for the feature weight updates.
p_init=0;                 %initial p.
p_max=0.1;                %maximum p.
p_step=0.01;              %p step.
t_max=100;                %maximum number of iterations.
beta_memory=0.2;          %amount of memory for the cluster weight updates.
fuzzy_degree=2;           %fuzzy membership degree

balance_parm = 0.1;       %balance parameter among to terms of loss function
T_pow = 1.4;              %power of T
a_coefficient = 3;        %coefficient of u
b_coefficient = 4;        %coefficient of t

I=0.00001;                %The value of this parameter is in the range of (0 and 1]
landa=I./var(X);
landa(landa==inf)=1;

window_size = 5;          % size of window for finding neighbors
NR = window_size^2 - 1;
beta_penalty = 0.8;
Neig = Find_Neighbors(window_size, reshape(X,size(f_ori,1),size(f_ori,2),d));

%Randomly initialize the cluster centers.
rng(1);
tmp=randperm(N);
M=X(tmp(1:k),:);

tic
[Cluster_elem,~,~,~]=mw_pfcm_brain_mri(X,M,k,p_init,p_max,p_step,t_max,beta_memory,N,fuzzy_degree,d,q,landa,balance_parm,a_coefficient,b_coefficient,T_pow, beta_penalty, Neig, NR);
Runtime = toc;

%% Lables (imagesize by imagesize)
[~,Lr2]=max(Cluster_elem);
Lr2 = reshape(Lr2,[size(f_ori,1) size(f_ori,2)]);
% result
a = double(reshape(GT, [1, size(GT,1) * size(GT,2)]));
b = double(reshape(Lr2, [1, size(GT,1) * size(GT,2)]));
EVAL = Evaluate(a',b');
accurcy=EVAL(1);
fm=EVAL(2);
nmi=EVAL(3);
Jaccard = EVAL(4);
Sensitivity = EVAL(5);
Lseg=Label_image(f_ori,Lr2);
figure,imshow(Lseg);
imwrite(Lseg, strcat('Results','\',FileName))
end