function new_RGB = myIncreaseSat(RGB,args)
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

% increase by % of max if vl is positive, and percentage of min if negative

if args.val>=0
    val = args.val*RGB(imax);
else
    val = args.val*RGB(imin);
end

% increase the Saturation 
if RGB(imax)==1
    %if max is already max possible increase only min 
    new_RGB(imin) = RGB(imin) + val;
else
    new_RGB(imax) = RGB(imax) + val;
    new_RGB(imin) = RGB(imin) + val;
end

% if any negative change to 0
new_RGB(new_RGB<0)=0;
new_RGB(new_RGB>1)=1;

end