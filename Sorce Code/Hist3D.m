function [ H ] = Hist3D( I )

H=zeros(256,256,256);
[U,V,Z] = size(I);

for u = 1:U
    for v = 1:V
        r = I(u,v,1);
        g = I(u,v,2);
        b = I(u,v,3);
        H(r+1,g+1,b+1) = H(r+1,g+1,b+1)+1;
    end
end

end

