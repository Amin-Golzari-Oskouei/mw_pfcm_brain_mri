clear ;close all;clc;

 path = 'Dataset';
imagefiles = dir([path,'\original','\*.tif']);

nfiles = length(imagefiles);    % Number of files found

m = {};
for ii = 1 : nfiles
    currentfilename = imagefiles(ii).name;
    currentfoldername = imagefiles(ii).folder;
    [accurcy,fm, nmi, Jaccard, Sensitivity] = main(currentfilename, path);
    [a,b] = max(fm);
    m{ii, 1} = currentfilename;
    m{ii, 2} = accurcy(b);
    m{ii, 3} = fm(b);
    m{ii, 4} = nmi(b);
    m{ii, 5} = Jaccard(b);
    m{ii, 6} = Sensitivity(b);
    writetable(table(m), 'Results.txt')
end
cellfun(@mean, m);     %provided you want a 2 x 82 numeric result.
sum(cellfun(@sum, m)) ./ sum(cellfun(@length, m))

