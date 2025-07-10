clear, close all, clc

path = 'Results\';

imagefiles = dir([path]);
nfiles = length(imagefiles);    % Number of files found
for ii = 3 : nfiles
    currentfilename = imagefiles(ii).name;
    currentfoldername = imagefiles(ii).folder;
    load([currentfoldername,'\',currentfilename]);
    currentfilename=currentfilename(1:end-4);
    mkdir('LabelsImg');
    IG = label2rgb(IG);
    imwrite(IG,['LabelsImg\',currentfilename ,'.bmp'])   
end