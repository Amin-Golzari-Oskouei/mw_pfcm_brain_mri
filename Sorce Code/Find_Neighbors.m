function result = Find_Neighbors(window_size, Img)
[row,col,~] =size(Img);
A = reshape([1:row*col],row,col);
% use the help of a bigger matrix
B=nan(size(A)+window_size-1);
B(ceil(window_size/2):end-floor(window_size/2),ceil(window_size/2):end-floor(window_size/2))=A;
B = fillmissing(B,'nearest');
B = fillmissing(B','nearest')';
% B(:,1:floor(window_size/2)) = B(:,ceil(window_size/2));
% B(:,end-floor(window_size/2:end))= B(:,end-ceil(window_size/2));
% pre-define memory for result
result = zeros(row*col,window_size*window_size);
% calculate!
n=1;
for i=ceil(window_size/2):size(A,1)+floor(window_size/2)
  for j=ceil(window_size/2):size(A,2)+floor(window_size/2)
    tmp=B(i-floor(window_size/2):i+floor(window_size/2),j-floor(window_size/2):j+floor(window_size/2));
    tmp(ceil(window_size/2),ceil(window_size/2))=nan;
    result(n,:)=reshape(tmp, 1,window_size*window_size);
    n=n+1;
  end
end

result(:,ceil(window_size*window_size/2)) = [];

end