function [lineOut, fillOut] = mySEMshade(amean,astd,alpha,acolor,F,smth,extra)
% usage: mySEMshade(M,S,alpha,acolor,F,smth)
% plot mean and sem/std shown as shading
% - M is the mean 
% - astd - is the SEM or STD (same size as M)
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% - method should be 'std' or 'sem' and defines the shading % added by magdalena 6/12/2023
%
% magdalena sabat 12/06/2024 adapted from smusall 2010/4/23


if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end

if exist('F','var')==0 || isempty(F)
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1; %no smoothing by default
end  

if ne(size(F,1),1)
    F=F';
end
amean = squeeze(amean);
astd = squeeze(astd);

if ne(size(amean,1),1)
    amean=amean';
end
if ne(size(astd,1),1)
    astd=astd';
end
% amean = nanmean(amatrix,1); %get man over first dimension
% if smth > 1
%     amean = boxFilter(nanmean(amatrix,1),smth); %use boxfilter to smooth data
% end

% switch method
%     case 'std'
%         astd = std(amatrix,[],1,'omitnan'); % to get std shading
%     case 'sem'
%         astd = std(amatrix,[],1,'omitnan')/sqrt(size(amatrix,1)); % to get sem shading
% end

if exist('extra')
    astd=mean(extra);
end

if exist('alpha','var')==0 || isempty(alpha) 
    fillOut = fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none');
    acolor='k';
else
    fillOut = fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none');
end

if ishold==0
    check=true; else check=false;
end

hold on;
lineOut = plot(F,amean,'Color',acolor,'linewidth',1); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end


function dataOut = boxFilter(dataIn, fWidth)
% apply 1-D boxcar filter for smoothing

fWidth = fWidth - 1 + mod(fWidth,2); %make sure filter length is odd
dataStart = cumsum(dataIn(1:fWidth-2),2);
dataStart = dataStart(1:2:end) ./ (1:2:(fWidth-2));
dataEnd = cumsum(dataIn(length(dataIn):-1:length(dataIn)-fWidth+3),2);
dataEnd = dataEnd(end:-2:1) ./ (fWidth-2:-2:1);
dataOut = conv(dataIn,ones(fWidth,1)/fWidth,'full');
dataOut = [dataStart,dataOut(fWidth:end-fWidth+1),dataEnd];

end

