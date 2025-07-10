function [Cluster_elem,M,W,Z]=mw_pfcm_brain_mri(X,M,k,p_init,p_max,p_step,t_max,beta_memory,N,fuzzy_degree,d,q,landa,balance_parm,a_coefficient,b_coefficient,T_pow, beta_penalty, Neig, NR)
%
%Arezou Najafi Moghaddam, Nasser Aghazadeh, Mahdi Hashemzadeh, Amin Golzari Oskouei, "An Adaptive Multi-Weighted Possibilistic Fuzzy C-Means
% Clustering Approach for Brain Tumor Segmentation in Magnetic Resonance Imaging", engineering of artificial intelligence, 2025.
%
%Function Inputs
%===============
%
%X is an Nxd data matrix, where each row corresponds to an instance.
%
%M is a kxd matrix of the initial cluster centers. Each row corresponds to a center.
%
%k is the number of clusters.
%
%p_init is the initial p value (0<=p_init<1).
%
%p_max is the maximum admissible p value (0<=p_max<1).
%
%p_step is the step for increasing p (p_step>=0).
%
%t_max is the maximum number of iterations.
%
%beta_memory controls the amount of memory for the weight updates (0<=beta<=1).
%
%fuzzy_degree is fuzzy membership degree
%
%balance_parm is balance parameter among to terms of loss function
%
%a_coefficient is coefficient of u
%
%b_coefficient is coefficient of b
% 
% T_pow is power of T
%
%Function Outputs
%================
%
%Cluster_elem is a kxd matrix containing the final cluster assignments.
%
%M is a kxd matrix of the final cluster centers. Each row corresponds to a center.
%
%W is a k-dimensional row vector containing the final weights of each cluster.
%
%z is a kxd matrix of the final weights of each fatuter in each cluster.
%
%Courtesy of A. Golzari Oskouei


if p_init<0 || p_init>=1
    error('p_init must take a value in [0,1)');
end

if p_max<0 || p_max>=1
    error('p_max must take a value in [0,1)');
end

if p_max<p_init
    error('p_max must be greater or equal to p_init');
end

if p_step<0
    error('p_step must be a non-negative number');
end

if beta_memory<0 || beta_memory>1
    error('beta must take a value in [0,1]');
end

if q==0
    error('beta must be a non-zero number');
end

if p_init==p_max
    
    if p_step~=0
        fprintf('p_step reset to zero, since p_max equals p_init\n\n');
    end
    
    p_flag=0;
    p_step=0;
    
elseif p_step==0
    
    if p_init~=p_max
        fprintf('p_max reset to equal p_init, since p_step=0\n\n');
    end
    
    p_flag=0;
    p_max=p_init;
    
else
    p_flag=1; %p_flag indicates whether p will be increased during the iterations.
end

%--------------------------------------------------------------------------
%Weights are uniformly initialized.
W=ones(1,k)/k;  %initial cluster weights
Z=ones(k,d)/d;  %initial faeture weights
T=ones(N,k)/k;  %initial faeture weights

%Other initializations.
p=p_init; %Initial p value.
p_prev=p-10^(-8); %Dummy value.
empty=0; %Count the number of iterations for which an empty or singleton cluster is detected.
Iter=1; %Number of iterations.
O_F_old=inf; %Previous iteration objective (used to check convergence).
Cluster_elem_history=[];
W_history=[];
z_history=[];
%--------------------------------------------------------------------------

fprintf('\nStart of fuzzy C-means clustering method based on feature-weight and cluster-weight learning iterations\n');
fprintf('----------------------------------\n\n');

%The proposed iterative procedure.
while 1
    %Update the cluster assignments.
    for j=1:k
        distance(j,:,:) = (1-exp((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2)));
        WBETA = transpose(Z(j,:).^q);
        WBETA(WBETA==inf)=0;
        dNK(:,j) = W(1,j).^p * reshape(distance(j,:,:),[N,d]) * WBETA   ;
        if Iter>1
            cc = (1-Cluster_elem(j,:)).^fuzzy_degree;
            dNK(:,j)= dNK(:,j) + (beta_penalty/NR).*sum(cc(Neig),2);
        end
        
    end
    
    tmp1 = zeros(N,k);
    for j=1:k
        tmp2 = (dNK./repmat(dNK(:,j),1,k)).^(1/(fuzzy_degree-1));
        tmp2(tmp2==inf)=0;
        tmp2(isnan(tmp2))=0;
        tmp1=tmp1+tmp2;
    end
    Cluster_elem = transpose(1./tmp1);
    
    
    Cluster_elem(isnan(Cluster_elem))=1;
    Cluster_elem(Cluster_elem==inf)=1;
    
    if nnz(dNK==0)>0
        for j=1:N
            if nnz(dNK(j,:)==0)>0
                Cluster_elem(find(dNK(j,:)==0),j) = 1/nnz(dNK(j,:)==0);
                Cluster_elem(find(dNK(j,:)~=0),j) = 0;
            end
        end
    end
    
    %Update gama.
    for j=1:k
        distance(j,:,:) = (1-exp((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2)));
        WBETA = transpose(Z(j,:).^q);
        WBETA(WBETA==inf)=0;
        dNK(:,j) = reshape(distance(j,:,:),[N,d]) * WBETA * W(1,j)^p ;
    end
    
    gama = balance_parm * 2 * sum(dNK .* transpose(Cluster_elem.^fuzzy_degree)) ./ sum(Cluster_elem.^fuzzy_degree,2)';
    
    %Update the T.
    T = (1+(repmat(2*b_coefficient./gama, N,1) .* dNK).^(1/(T_pow-1))).^-1;
    
    %Calculate the fuzzy C-means clustering method based on feature-weight and cluster-weight learning objective.
    O_F=object_fun(N,d,k,Cluster_elem,landa,M,fuzzy_degree,W,Z,q,p,X,gama,T,a_coefficient,b_coefficient,T_pow,beta_penalty, Neig, NR);
    
    %If empty or singleton clusters are detected after the update.
    for i=1:k
        
        I=find(Cluster_elem(i,:)<=0.05);
        if length(I)==N-1 || length(I)==N
            
            fprintf('Empty or singleton clusters detected for p=%g.\n',p);
            fprintf('Reverting to previous p value.\n\n');
            
            O_F=NaN; %Current objective undefined.
            empty=empty+1;
            
            %Reduce p when empty or singleton clusters are detected.
            if empty>1
                p=p-p_step;
                
                %The last p increase may not correspond to a complete p_step,
                %if the difference p_max-p_init is not an exact multiple of p_step.
            else
                p=p_prev;
            end
            
            p_flag=0; %Never increase p again.
            
            %p is not allowed to drop out of the given range.
            if p<p_init || p_step==0
                
                fprintf('\n+++++++++++++++++++++++++++++++++++++++++\n\n');
                fprintf('p cannot be decreased further.\n');
                fprintf('Either p_step=0 or p_init already reached.\n');
                fprintf('Aborting Execution\n');
                fprintf('\n+++++++++++++++++++++++++++++++++++++++++\n\n');
                
                %Return NaN to indicate that no solution is produced.
                M=NaN(k,size(X,2));
                return;
            end
            
            %Continue from the assignments and the weights corresponding
            %to the decreased p value.
            a=(k*empty)-(k-1);
            b=k*empty;
            Cluster_elem=Cluster_elem_history(a:b,:);
            W=W_history(empty,:);
            aa=(k*empty)-(k-1);
            bb=k*empty;
            Z=z_history(aa:bb,:);
            break;
        end
    end
    
    if ~isnan(O_F)
        fprintf('p=%g\n',p);
        fprintf('The clustering objective function is %f\n\n',O_F);
    end
    
    %Check for convergence. Never converge if in the current (or previous)
    %iteration empty or singleton clusters were detected.
    if Iter>=t_max || ~isnan(O_F) && ~isnan(O_F_old) && (abs(1-O_F/O_F_old) < 1e-6 )
        
        fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
        fprintf('Converging for p=%g after %d iterations.\n',p,Iter);
        fprintf('The final objective function is =%f.\n',O_F);
        fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
        
        break;
        
    end
    
    O_F_old=O_F;
    
    %Update the cluster centers.
    mf = Cluster_elem.^fuzzy_degree;       % MF matrix after exponential modification
    tf = T'.^T_pow;
    
    for j=1:k
        M(j,:) = (((a_coefficient.*(mf(j,:))) + (b_coefficient.*(tf(j,:)))) * (X .* (exp((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2)))))./...
            (((((a_coefficient.*(mf(j,:))) + (b_coefficient.*(tf(j,:))))*(exp((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2)))))); %new center
    end
    
    %Increase the p value.
    if p_flag==1
        
        %Keep the assignments-weights corresponding to the current p.
        %These are needed when empty or singleton clusters are found in
        %subsequent iterations.
        Cluster_elem_history=[Cluster_elem;Cluster_elem_history];
        W_history=[W;W_history];
        z_history=[Z;z_history];
        
        p_prev=p;
        p=p+p_step;
        
        if p>=p_max
            p=p_max;
            p_flag=0;
            fprintf('p_max reached\n\n');
        end
    end
    
    W_old=W;
    z_old=Z;
    
    %Update the feature weights.
    
    for j=1:k
        distance(j,:,:) = (1-exp((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2)));
        dWkm(j,:) = ((a_coefficient.*(Cluster_elem(j,:).^fuzzy_degree)) + (b_coefficient.*(T(:,j).^T_pow))') * reshape(distance(j,:,:),[N,d]);
    end
    
    tmp1 = zeros(k,d);
    for j=1:d
        tmp2 = (dWkm./repmat(dWkm(:,j),1,d)).^(1/(q-1));
        tmp2(tmp2==inf)=0;
        tmp2(isnan(tmp2))=0;
        tmp1=tmp1+tmp2;
    end
    Z = 1./tmp1;
    Z(isnan(Z))=1;
    Z(Z==inf)=1;
    
    if nnz(dWkm==0)>0
        for j=1:k
            if nnz(dWkm(j,:)==0)>0
                Z(j,find(dWkm(j,:)==0)) = 1/nnz(dWkm(j,:)==0);
                Z(j,find(dWkm(j,:)~=0)) = 0;
            end
        end
    end
    
    %Update the cluster weights.
    for j=1:k
        distance(j,:,:) = (1-exp((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2)));
        WBETA = transpose(Z(j,:).^q);
        WBETA(WBETA==inf)=0;
        Dw(1,j) = transpose(WBETA) * transpose(reshape(distance(j,:,:),[N,d])) *  transpose(((a_coefficient.*(Cluster_elem(j,:).^fuzzy_degree)) + (b_coefficient.*(T(:,j).^T_pow))')) ;
    end
    
    tmp = sum((repmat(Dw,k,1)./transpose(repmat(Dw,k,1))).^(1/(p-1)));
    tmp(tmp==inf)=0;
    tmp(isnan(tmp))=0;
    W = 1./tmp;
    W(isnan(W))=1;
    W(W==inf)=1;
    
    if nnz(Dw==0)>0
        W(find(Dw==0)) = 1/nnz(Dw==0);
        W(find(Dw~=0)) = 0;
    end
    
    %Memory effect.
    W=(1-beta_memory)*W+beta_memory*W_old;
    Z=(1-beta_memory)*Z+beta_memory*z_old;
    
    Iter=Iter+1;
end
end



