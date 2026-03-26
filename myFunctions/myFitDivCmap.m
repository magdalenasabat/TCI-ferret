function cmap = myFitDivCmap(data, colors)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates colormap centered around zero
%
% date: 06/10/2023
% author: Magdalena Sabat magdalena.sabat@ens.psl.eu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('colors','var')
  colors =[0    0   .5; ...
           1,   1,  1; ...
           .5   0   0];
end

% define center
if any(data > 0,'all') && any(data < 0,'all')
    center = 0;
else
    warning('No pos&&neg values in data. Returning Cmap centered around mean value')
    center = mean(data,'all'); 
end

% define relative distance to center
dist_neg = abs(min(data,[],'all')) / (abs(min(data,[],'all')) + abs(max(data,[],'all')));

% create colormap
n = floor(256*dist_neg);

half1 = [   [linspace(colors(1,1),colors(2,1),n) ]; ...
            [linspace(colors(1,2),colors(2,2),n) ]; ...
            [linspace(colors(1,3),colors(2,3),n) ] ];
            
half2 = [   [linspace(colors(2,1),colors(3,1),256 - n) ]; ...
            [linspace(colors(2,2),colors(3,2),256 - n) ]; ...
            [linspace(colors(2,3),colors(3,3),256 - n) ] ];

cmap = [half1,half2]';

end