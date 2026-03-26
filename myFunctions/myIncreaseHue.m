function new_RGB = myIncreaseHue(RGB,args)
% increases hue of the RGB triplet
% magdalena sabat 25/2024 magdalenajoannasbat@gmail.com
arguments
    RGB (1,3) double
    args.val (1,1) double = .1 % increase in hue, in precentage of max
end

%To saturate a color we need to increase the difference between the lowest and highest RGB value
[~,imax] = max(RGB);
[~,imin] = min(RGB);

new_RGB=RGB;

% increase the Saturation 
new_RGB(imax) = RGB(imax) + args.val*RGB(imax);
new_RGB(imin) = RGB(imin) + args.val*RGB(imax);
end