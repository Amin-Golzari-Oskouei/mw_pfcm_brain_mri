function [label,rgb,Eticets] = Coloring(I, C)

% input:
%         img: N * d where N is the number of pixels & d is the number of color channels
%         M:  k, d where k is the number of clusters and d is the number of color channels

% output:
%         distance: K * N where k is the number of clusters & N is the number of pixels

[u v z] = size(I);
img = double(reshape(I,u*v,3))/255;
M = C/255;
for j=1:size(M, 1)
    center2mat = double(repmat(M(j,:),size(img, 1),1));
    distance(j,:) = sum(1-exp((-1).*((img-center2mat).^2)),2);
    %     distance(j,:) = sqrt(sum((img-center2mat).^2,2));
end

[~, ind] = min(distance);
IG = M(ind,:);
label = reshape(ind,u,v,1);
rgb = reshape(IG,u,v,z);
Eticets = unique(ind);

% R = (I(:,:,1)); G = int32(I(:,:,2)); B = int32(I(:,:,3));
% center = C';
% MeanR =  center(1,:)';
% MeanG =  center(2,:)';
% MeanB =  center(3,:)';
% for i = 1:numel(Eticets)
%     newCentR(i) = round(mean(R(label==i)));
%     newCentG(i) = round(mean(G(label==i)));
%     newCentB(i) = round(mean(B(label==i)));
% end
% 
% newC = [newCentR;newCentG;newCentB]';


end



