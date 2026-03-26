function draw_connectedlines_violin(violindata, plotcolor, plottags)
% T Engelen March 2024

% function for plotting paired data in 2 conditions, for each showing the
% distribution and connected by lines to see change between conditions.

% Input:
% violindata: structure where each column represents inputs of a condition,
% each row is a participant (! function specifically made for 2 conditions)
% plotcolor: specified in RGB triplet
% plottags, two strings representing the x tags in the plot

nviolins = size(violindata,2); % should be two

for iviolin = 1:nviolins
    datavect  = violindata(:,iviolin);
    if size(plotcolor,1)
        colorvect = plotcolor(1,:);
    else
        colorvect = plotcolor(iviolin,:);
    end
    position  = iviolin;
 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% (1) draw violin shape:
    % computing the kernel density for drawing the violin
    if min(datavect) == max(datavect) % check if data is not distributed
        density = repelem(1, size(datavect,1));
        value   = datavect';
    elseif min(datavect) ~= max(datavect)
        if std(datavect) > realmin && iqr(datavect) > realmin
            [density, value] = ksdensity(datavect, 'bandwidth', 0.9 * min(std(datavect), iqr(datavect)/1.34) * numel(datavect)^(-1/5) );
            density          = density(value >= min(datavect) & value <= max(datavect));
        elseif std(datavect) < realmin
            [density, value] = ksdensity(datavect, 'bandwidth', 0.9 * (iqr(datavect)/1.34) * numel(datavect)^(-1/5) );
            density          = density(value >= min(datavect) & value <= max(datavect));
        elseif iqr(datavect) < realmin
            [density, value] = ksdensity(datavect, 'bandwidth', 0.9 * std(datavect) * numel(datavect)^(-1/5) );
            density          = density(value >= min(datavect) & value <= max(datavect));
        end
    end
    

    if sum(density) == 0
        warning('bandwidth formula fails => setting arbitrarily => check distribution')
        [density, value] = ksdensity(datavect, 'bandwidth', 0.05);
        density          = density(value >= min(datavect) & value <= max(datavect));
    end
    

    value      = value(value >= min(datavect) & value <= max(datavect));
    value(1)   = min(datavect);
    value(end) = max(datavect);
    width      = 0.3/max(density);
    
    % draw violin shape using patch
    hold on
    if iviolin == 1 % violin 1 is plotted pointing to the left
        X_patch = [repelem(position-0.05, size(density,2)) (position-0.05)-density(end:-1:1)*width];
    else  % violin 2 is plotted pointing to the right
        X_patch = [repelem(position+0.05, size(density,2)) (position+0.05)+density(end:-1:1)*width];

    end

    Y_patch = [value value(end:-1:1)];
    violins(position).violinshape = fill(X_patch, Y_patch , 1, 'facecolor', 0.5*(colorvect+1), 'edgecolor','none', 'facealpha', 0.5, 'LineWidth', 0.5);
    
    
end

% plot connected lines with alpha value
alphacolor = plotcolor;
alphacolor(4) = 0.6; % alpha value
plot(violindata', '-','Color',alphacolor,'LineWidth', 1.2); 
hold on
% using scatter to plot individual dots so we can set alpha values
scatter(ones(size(violindata,1),1),violindata(:,1)', 65,'filled', 'MarkerFaceColor', plotcolor, 'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75)
scatter(ones(size(violindata,1),1)+1,violindata(:,2)', 65,'filled', 'MarkerFaceColor', plotcolor, 'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75)


set(gca, 'box', 'off', 'TickDir','out', 'FontSize', 12, 'XTick', 1:nviolins, 'XTickLabel', plottags, 'XLim', [0.5 2.5],...
    'PlotBoxAspectRatio',[1,1,1], 'linewidth', 1);


