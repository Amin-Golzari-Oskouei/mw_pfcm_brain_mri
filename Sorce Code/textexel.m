clear close all clc
a = xlsread('d:\\3dHist.xlsx');
b = a(:,1);
c = a(:,2:end);
isn = ~isnan(b(:,1));
[u v] = find(isn);
p = zeros(600,15);
iter = 1;
for i = 1:length(u)-1
%     p = a(u(i,1));
    [x y] = size(c(u(i):u(i+1)-1,:)');
    p(iter,1) = a(u(i,1));
    p(iter:iter+2,2:y+1) =  c(u(i):u(i+1)-1,:)';
    iter = iter +3;
end
i = i+1;
[x y] = size(c(u(i):end,:)');
p(iter,1) = a(u(i,1));
p(iter:iter+2,2:y+1) =  c(u(i):end,:)';

xlswrite('d:\\Book2.xlsx',p);