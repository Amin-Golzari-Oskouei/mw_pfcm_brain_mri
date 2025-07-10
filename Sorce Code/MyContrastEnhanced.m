function [D] = MyContrastEnhanced(D,no)
%% initialization
[n1,n2]=size(D);
M=std(D(:));
p=0.5;
d=0;
k=3;
m=k*M;
%% Ranges and Conditions
for k1=1:n1
    for k2=1:n2
        if D(k1,k2)<no
            s=1;
            D(k1,k2)=s*D(k1,k2);
        elseif (no<D(k1,k2))&&(D(k1,k2)<2*no)
            s=(((D(k1,k2)-no)/no)*((m/no)^p))+(((2*no)-D(k1,k2))/no);
            D(k1,k2)=s*abs(D(k1,k2));
        elseif (2*no<D(k1,k2))&&(D(k1,k2)<m)
            s=(m/D(k1,k2))^p;
            D(k1,k2)=s*abs(D(k1,k2));
        elseif m<=(D(k1,k2))
            s=(m/D(k1,k2))^d;
            D(k1,k2)=s*abs(D(k1,k2));
        end
    end
end
end