function [smoothedRGB] = NonLocalMeansFilter(f_ori)
    % تبدیل تصویر RGB به خاکستری برای برآورد نویز
    grayImage = rgb2gray(f_ori);

    % برآورد نویز برای تنظیم پارامترهای NLM
    estimatedSigma = std2(double(grayImage));

    % اعمال فیلتر NLM به کانال‌های رنگ به صورت جداگانه
    smoothedR = imnlmfilt(f_ori(:,:,1), 'DegreeOfSmoothing', 10*estimatedSigma^2);
    smoothedG = imnlmfilt(f_ori(:,:,2), 'DegreeOfSmoothing', 10*estimatedSigma^2);
    smoothedB = imnlmfilt(f_ori(:,:,3), 'DegreeOfSmoothing', 10*estimatedSigma^2);

    % ترکیب کانال‌ها
    smoothedRGB = cat(3, smoothedR, smoothedG, smoothedB);
end
