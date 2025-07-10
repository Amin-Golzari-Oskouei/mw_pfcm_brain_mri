function [Output_Image] = EdgeEnhancedCurvelet(f_ori)
%Read Input Image
Input_Image=f_ori;
%Red Component of Colour Image
Red_Input_Image=Input_Image(:,:,1);
%Green Component of Colour Image
Green_Input_Image=Input_Image(:,:,2);
%Blue Component of Colour Image
Blue_Input_Image=Input_Image(:,:,3);

%% Apply Curvelet On Red chanel
X_Red=double(Red_Input_Image);
[m1,n1]=size(X_Red);
%take curvelet transform
C_Red = fdct_wrapping(X_Red,1,1);
Ct_Red=C_Red;
% coarse approximation is set to zero
for r_Red=1:length(C_Red{end})
    C_Red{end}{r_Red}=zeros(size(C_Red{end}{r_Red}));
end
% multiply other subbands by amplification factor
for w_Red=1:(length(C_Red)-1)
    for s_Red=1:length(C_Red{w_Red})
        C_Red{w_Red}{s_Red}=2*C_Red{w_Red}{s_Red};
    end
end
%inverse curvelet transform
First_Level_Decomposition(:,:,1) = real(ifdct_wrapping(C_Red,1,m1,n1));  %Y is changed img
% Output_Image=uint8(First_Level_Decomposition(:,:,1));
% subplotting results
%    subplot(1,2,1);imshow(Input_Image,[]);title('Original image');
%    subplot(1,2,2);imshow(Output_Image,[]);title('Enhanced image');
%% Apply Curvelet On Green chanel
X_Green=double(Green_Input_Image);
[m2,n2]=size(X_Green);
C_Green = fdct_wrapping(X_Green,1,1);
Ct_Green=C_Green;
% coarse approximation is set to zero
for r_Green=1:length(C_Green{end})
    C_Green{end}{r_Green}=zeros(size(C_Green{end}{r_Green}));
end
% multiply other subbands by amplification factor
for w_Green=1:(length(C_Green)-1)
    for s_Green=1:length(C_Green{w_Green})
        C_Green{w_Green}{s_Green}=2*C_Green{w_Green}{s_Green};
    end
end
%inverse curvelet transform
 First_Level_Decomposition(:,:,2) = real(ifdct_wrapping(C_Green,1,m2,n2));
 %% Apply Curvelet On Blue chanel

 X_Blue=double(Blue_Input_Image);
 [m3,n3]=size(X_Blue);
 C_Blue = fdct_wrapping(X_Blue,1,1);
 Ct_Blue=C_Blue;
% coarse approximation is set to zero
for r_Blue=1:length(C_Blue{end})
    C_Blue{end}{r_Blue}=zeros(size(C_Blue{end}{r_Blue}));
end
% multiply other subbands by amplification factor
for w_Blue=1:(length(C_Blue)-1)
    for s_Blue=1:length(C_Blue{w_Blue})
        C_Blue{w_Blue}{s_Blue}=2*C_Blue{w_Blue}{s_Blue};
    end
end
%inverse curvelet transform
 First_Level_Decomposition(:,:,3) = real(ifdct_wrapping(C_Blue,1,m3,n3));
 Output_Image=uint8(First_Level_Decomposition);

end