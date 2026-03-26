% Loads published data and reproduces correlation computation and model
% fitting.
% 
% Sabat, M., Gouyette, H., Gaucher, Q., Espejo, M. L., David, S. V., Norman-Haignere, S., & Boubenec, Y. (2025). 
% Neurons in auditory cortex integrate information within a constrained and context-invariant temporal window. 
% Current Biology, 0(0). https://doi.org/10.1016/j.cub.2025.11.011
% 
%
% author: Magdalena email: magdalena.sabat@ens.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize
clc; close all; clear all; 


addpath(genpath(strcat([ 'code' filesep 'myFunctions'])));
addpath(genpath(strcat([ 'code' filesep 'toolbox-sam-tci'])));

dataFolder = strcat([ 'data' filesep 'experiment1']);


subnames={'A','B','C'};
sessions = {[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,26,27,29,30,31], ...
    [1,2,3,4,5,7,8,9,10,12,13,14,15,22,23,25,26,27,28,30,31,32,34,35,38],...
    [2,3,4,5,6,7,9,10,11,12,13,14,15]};

% for plotting 
markerstyles = {'^','square', 'o'}; 

%% Example for the 1st session of animal A

xi=1;
all_widths = zeros(1,500);
all_distances = zeros(1,500);
all_subs = zeros(1,500);
all_areas = cell(1,500);
all_CCC = zeros(696,10,500);
all_NC = zeros(696,10,500);
 


for subid= 1:length(subnames) 
    
    subjid = subnames{subid};
    fprintf(strcat([' SUBJECT ' subjid '\n']))
    
    
    for sessN=sessions{subid}
        
%         fprintf(strcat([' sess #' sprintf('%03d',sessN) '\n']))

        %load files
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_anatomical.mat'  ]));
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_model_fits.mat'  ]));
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_cross_context_correlation.mat'  ]));
 

        % pickup model fits and distance
        units_n = length(M.best_intper_sec);
        all_widths(xi:xi+units_n-1)=M.best_delay_sec_median;
        all_distances(xi:xi+units_n-1)=distances;
        all_subs(xi:xi+units_n-1)=subid;
        all_areas(xi:xi+units_n-1)=areas;
        
        all_CCC(:,:,xi:xi+units_n-1)=L.diff_context;
        all_NC(:,:,xi:xi+units_n-1)=L.same_context;
        xi=xi+units_n;
        
    end
end


all_widths = all_widths(1:xi-1); all_widths=all_widths*1000;
all_distances = all_distances(1:xi-1);
all_areas = all_areas(1:xi-1);
all_subs = all_subs(1:xi-1);
all_CCC = all_CCC(:,:,1:xi-1);
all_NC = all_NC(:,:,1:xi-1);
lag_t = L.lag_t;
unique_segs = L.unique_segs;
%% Fig 1c


f = figure('Units','centimeters','Position',[100 100 5 5]); 
ax = gca;
axis(ax, 'square'); hold on

for subid= 1:length(subnames) 

    markerstyle = markerstyles{subid};
    mask = all_subs==subid;
    
    scatter(ax,all_distances(mask),all_widths(mask),5,'marker',markerstyle, ...
                'MarkerEdgeColor',[.2 .2 .2],...
                'MarkerEdgeAlpha' ,.7,...
                'MarkerFaceAlpha',.3,...
                'jitter','on','JitterAmount',.3,'HandleVisibility', 'off');
end

[all_dist_sorted,ind_dist] = sort(all_distances);

[W,Sw] =   polyfit( all_dist_sorted, log2(all_widths(ind_dist)),1);
[ylw, deltaw] = polyval(W,all_dist_sorted,Sw);


plot(ax,all_dist_sorted,2.^ylw,'--','LineWidth',2, 'Color',[0.8500 0.3250 0.0980] ,'DisplayName', 'regression');
text( ax, median(all_dist_sorted)+.15,median(2.^ylw)-23, [ '$\rho = ' num2str(round(rankcorr(all_distances',all_widths'),2)) '$'], 'Color',[0.8500 0.3250 0.0980],'Interpreter', 'latex')

set(ax, 'YScale', 'log','YMinorTick','off','YTick',[ 8 16 32 64 128],...
    'ylim', [10 150], 'xlim', [0 6], 'XTick',[0 2 4 6] )

ylabel(ax, 'width (ms)')
xlabel(ax, 'distance to PAC (mm)')

set(findall(ax,'-property','FontSize'),'FontSize',11)
%% Fig 1d


hist_col = [[64 43 183];[51 107 51]]/255;

figure('Units', 'centimeters','Position', [50,50,3, 2.65]);
hold on;
bins = 10.^(1:.05:2.4);

mov=[42 35];
mov=[28 24];
ticks = [8  32  128 ];
gr=unique(all_areas);

for j=1:length(gr)
    area_mask = strcmp(all_areas, gr(j));
    d = all_widths(area_mask);

    h=histogram(d,bins,'FaceColor',hist_col(j,:),'EdgeColor',hist_col(j,:),'FaceAlpha',.7); %',.4,[.5 .5 .5],lag_t(2:end-1),[], 'sem');
    
    xlim([min(ticks) max(ticks)+30])
    set(gca,'XScale','log')%
    xtickangle(0)
    if j ==2
        set(gca, 'xtick', ticks, 'XTickLabels', ticks,'XMinorTick','off')
    else
        set(gca, 'xtick', ticks, 'XTickLabels', [],'XMinorTick','off')
    end
    set(gca, 'ytick', [0:10:100],'Ylim',[0 max(h.Values)*1.5])


    xline(median(d),'Color',hist_col(j,:))
    text(median(d)-11,mov(j),strcat(['Median=' num2str(round(median(d),1))]),'FontSize',6 ,'Color',hist_col(j,:))   
   
    
    
end
xlabel('width')  
ylabel('Neurons/bin')
    
set(findall(gcf,'-property','FontSize'),'FontSize',7)

%% 

hist_col = [[64 43 183];[51 107 51]]/255;
valid_durs = [1,3,4,5,7,9];
axgrid = [1,numel(valid_durs)];  % [#rows, #cols]
titles = {'closest', 'intermediate', 'farthest'}; 

figTCI = figure('Units', 'centimeters','Position', [50,50,9.9,2.7]);
% tclMain = tiledlayout(axgrid(1),1); 
% tcl = gobjects(1,axgrid(1));
% ax = gobjects(axgrid); 

titles={'Primary Auditory Cortex - MEG','Non-primary Auditory Cortex - PEG'};
idexes = [1:1:axgrid(2)];

axs = tiledlayout(axgrid(1),axgrid(2))

k=1;
gr = unique(all_areas);
for i= 1:axgrid(2)
%     if i==axgrid(2); continue; end
    idx = valid_durs(i);

    ax(k,i) = nexttile(axgrid(2)*(k-1)+i);
    hold( ax(k,i) , 'on')
    for j=1:length(gr)
        dists_mask = strcmp(all_areas, gr(j));
        s = squeeze(all_NC(:,idx,dists_mask));
        d =squeeze( all_CCC(:,idx,dists_mask));


    %         semshade(s(2:end-1,dists_mask)',.4,[.5 .5 .5],lag_t(2:end-1),[], 'std');
    %         plot(lag_t(2:end-1),median([d([2:end-1],dists_mask)./s(2:end-1,dists_mask)]'),'color',hist_col(j,:));
        plot(lag_t,median(d./s,2),'LineWidth',0.5,'color',hist_col(j,:))
        yline(1,'color',[.5 .5 .5])

        xlim([-.05 .25]);
        ylim([-0.15 1.15])

        ix=[0 unique_segs(idx)/1000 unique_segs(idx)/1000 0 ];
        yl = ylim;
        iy = [yl(1) yl(1) yl(2) yl(2)];
        p = patch('XData',ix,'YData',iy, 'FaceColor',[.85 .85 .85], 'LineStyle', 'none');hold on;
        uistack(p,'bottom')

        if i==5 & j==3; xlabel('time(ms)'); end

        max_idx = find(lag_t(2:end-1)>=unique_segs(idx)/1000,1,'first');
%         dist_same = s(max_idx,dists_mask)';
%         dist_diff = d(max_idx,dists_mask)';
    %         overlap = normalized_overlap(dist_same,dist_diff)

        if i==1
            yticks([0 1])
            title({'     '})
    %             yticklabels({'0' 'max'})

        else
            yticks([])
            yticklabels([])  
        end

        if j==2 & i==1
            xticks([0  .2])
            xtickangle(0)
            xticklabels([0 200])
        else
            xticks([0 .2])
            xtickangle(0)
            xticklabels([])
        end
    end

end

set(findall(figTCI,'-property','FontSize'),'FontSize',7)