function [ F,F2,Q ] = Performance_Eval( m,I2,Label,center )
%PERFORMANCE_EVAL Summary of this function goes here
%   Detailed explanation goes here
disp('Performance evaluation ...')

global I
tttt = [nonzeros(unique(Label))'];
y = 0;
R = int32(I(:,:,1)); G = int32(I(:,:,2)); B = int32(I(:,:,3));
MeanR = []; MeanG = []; MeanB = [];
AreaSeg = []; 
EuclidnDistR = []; EuclidnDistG = []; EuclidnDistB = [];
RA = [];
MeanR =  center(1,:)';
MeanG =  center(2,:)';
MeanB =  center(3,:)';
for i = tttt
    y = y+1;  
    AreaSeg = [AreaSeg, numel(Label(Label==i))]; 
%     EuclidnDistR =  (R(Label==i) - round(MeanR(i))).^2;
%     EuclidnDistG =  (G(Label==i) - round(MeanG(i))).^2;
%     EuclidnDistB =  (B(Label==i) - round(MeanB(i))).^2;
        newCentR(i) = round(mean(R(Label==i)));
        newCentG(i) = round(mean(G(Label==i)));
        newCentB(i) = round(mean(B(Label==i)));
        EuclidnDistR =  (R(Label==i) - newCentR(i)).^2;
        EuclidnDistG =  (G(Label==i) - newCentG(i)).^2;
        EuclidnDistB =  (B(Label==i) - newCentB(i)).^2;
    EuclidnDist(y) = sum(EuclidnDistR);        
    EuclidnDist(y) = sum(sqrt(double(EuclidnDistR + EuclidnDistG + EuclidnDistB)))/AreaSeg(i);   
    RA(y) = size(AreaSeg(AreaSeg==AreaSeg(y)),2);
    
end

% [w,h] = size(I(:,:,1));
[U,V,Z] = size(I);
F = (1/(U*V*1000)) * sqrt(m) * sum((EuclidnDist.^2)./sqrt(double(AreaSeg)));
Q = (1/(U*V*10000)) * sqrt(m) * sum((EuclidnDist.^2)./((1+log10(AreaSeg)))+ ((RA)./(AreaSeg)).^2);
tepm = [];
for i = 1:m
    temp(i) = RA(i)^(1+(1/AreaSeg(i)));
end
F2 = (1/(U*V*10000)) * sqrt(sum(temp))*sum((EuclidnDist.^2)./sqrt(double(AreaSeg)));

disp(['F = ',sprintf('%.8f ', F)])
disp(['F" = ',sprintf('%.8f ', F2)])
disp(['Q = ',sprintf('%.8f ', Q)])
end

