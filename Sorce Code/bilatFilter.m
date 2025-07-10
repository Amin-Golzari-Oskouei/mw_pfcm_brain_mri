function [smoothedRBG2] = bilatFilter(f_ori)
imRGB = f_ori;
% figure,imshow(imRGB);
imLAB = rgb2lab(imRGB);
patch = imLAB;
patchSq = patch.^2;
edist = sqrt(sum(patchSq,3));
patchVar = std2(edist).^2;
DoS2 = 0.9*patchVar;
sigma = 4;
smoothedLAB2 = imbilatfilt(imLAB,DoS2,sigma);
smoothedRBG2 = lab2rgb(smoothedLAB2,"Out","uint8");
end