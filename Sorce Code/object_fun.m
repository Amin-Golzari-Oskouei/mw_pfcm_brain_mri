function j_fun = object_fun(N,d,k,Cluster_elem,landa,M,fuzzy_degree,W,z,beta_z,p,X,gama,T,a_coefficient,b_coefficient,T_pow,beta_penalty, Neig, NR)

term_3 = 0;
for j=1:k
    distance(j,:,:) = (1-exp((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2)));
    WBETA = transpose(z(j,:).^beta_z);
    WBETA(WBETA==inf)=0;
    dNK(:,j) = reshape(distance(j,:,:),[N,d]) * WBETA * W(1,j)^p ;
    
    %term 3
    cc = (1-Cluster_elem(j,:)).^fuzzy_degree;
    term_3=term_3+sum((Cluster_elem(j,:).^fuzzy_degree)'.* sum(cc(Neig),2));
    
end
term_1 = 2 * sum(sum(dNK .* ((a_coefficient.*transpose(Cluster_elem.^fuzzy_degree)) + (b_coefficient.*(T.^T_pow)))));
term_2 = sum(sum((1-T).^T_pow).* gama);
j_fun = term_1+term_2+(beta_penalty/NR)*term_3;
end

