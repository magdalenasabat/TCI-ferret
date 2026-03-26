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


addpath(genpath(strcat(['code' filesep 'myFunctions'])));
addpath(genpath(strcat(['code' filesep 'toolbox-sam-tci'])));

dataFolder = strcat([ 'data' filesep 'experiment2']);

subnames={'D','E','F','G','H'};
sessions = {[17:29,31,33], ...
    [29:32],...
    [3:5,11,12,18,19],...
    [9,10,15:18,23,24],...
    [10,21,23]};

% for plotting 
markerstyles = {'^','square', 'o'}; 

%% Get all data

xi=1;
all_widths = zeros(1,500);
all_layers = zeros(1,500);
all_areas = cell(1,500);
all_subs = zeros(1,500);
all_CCC = zeros(691,6,500);
all_NC = zeros(691,6,500);
 


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
        all_widths(xi:xi+units_n-1)=M.best_intper_sec;
        all_subs(xi:xi+units_n-1)=subid;
        all_layers(xi:xi+units_n-1)=layer;
        all_areas(xi:xi+units_n-1)=area;
        
        all_CCC(:,:,xi:xi+units_n-1)=L.diff_context;
        all_NC(:,:,xi:xi+units_n-1)=L.same_context;
        xi=xi+units_n;
        
    end
end


all_widths = all_widths(1:xi-1); all_widths=all_widths*1000;
all_layers = all_layers(1:xi-1);
all_areas = all_areas(1:xi-1);
all_subs = all_subs(1:xi-1);
all_CCC = all_CCC(:,:,1:xi-1);
all_NC = all_NC(:,:,1:xi-1);
lag_t = L.lag_t;
unique_segs = L.unique_segs;

%% Fig 2b cumulative distributiotns

ticks = [ 32  128  512];
cols = lines(3);
% cols=[cols(1,:),cols(1,:),cols(1,:)]
gr=unique(all_areas);
lr = [4 13 56]; % layers: 13-> supra 56-> infra 4-> granular
for ar =gr

   
    figure('Units','centimeters','Position',[0 0 4 6  ]); 
    ax1 = axes(gcf);  hold on
    axis(ax1, 'square'); 
    hold on;

    si=1;
    s=1;
    for sx=lr
    
        mask = strcmp(all_areas,ar ) & (all_layers==sx) ;

          % cumulative distribbution
        [f,x,flo,fup] = ecdf(all_widths(mask),'Bounds','on','Alpha',0.01);
        flo=fillmissing(flo,'linear');
        fup=fillmissing(fup,'linear');
   
        plot(x,f,'LineWidth',.5,'color',cols(s,:));
        fill([x' fliplr(x')  ],[flo' fliplr(fup') ],cols(s,:),'FaceAlpha',.1,'linestyle','none',HandleVisibility='off')
      
        
        set(gca,'XScale','log')
%         ticks = [7 15 32 64 128 256 512];
        xlim([min(ticks)-18 max(ticks)]);
        ylim([0 1])
        set(gca, 'xtick', ticks, 'XTickLabels', ticks,'XMinorTick','off',...
            'ytick', [0 .5 1], 'YTickLabels',[])
        s=s+1;
        
    end
    xtickangle(0)
    xlabel('width')
    si=si+1;
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    title(ar)
    legend('granular','supragranular','infragranular','Location','southoutside','box','off')
%     export_fig( [tt filesep 'cumulative_width_' ar{:} '.png'],gcf,'-r2000' )
end


%% Fig 2c CCC


valid_durs = [1:6];
axgrid = [1,numel(valid_durs)];  % [#rows, #cols]

titles = {'A1', 'PEG'}; 
for li=1:length(titles)
    
    figTCI = figure('Units', 'centimeters','Position', [0 0 9 3]);
    axs = tiledlayout(axgrid(1),axgrid(2),'TileSpacing' , 'compact');
    sgtitle(titles(li));
    idexes = [1:1:axgrid(2)];
    ar=gr(li);
    
    for i= 1:axgrid(2)
   
        
        
        ax(i) =  nexttile(i);
         for j=1:numel(lr)
             sx=lr(j);
            idx = valid_durs(i);
            mask = strcmp(all_areas,ar ) & (all_layers==sx) ;

            hold( ax(i) , 'on')
          
            x1=median(squeeze(all_NC(:,i,mask)),2);
            x2=median(squeeze(all_CCC(:,i,mask)),2);
            plot(lag_t,x2./x1,'color',cols(j,:));   
            yline(1)

            xlim([-.05 .25]);
            ylim([-0.25 1.5])
            title(' ')
            ax(i).TitleFontSizeMultiplier = .2;
            ix=[0 unique_segs(idx)/1000 unique_segs(idx)/1000 0 ];
            yl = ylim;
            iy = [yl(1) yl(1) yl(2) yl(2)];
            p = patch('XData',ix,'YData',iy, 'FaceColor',[.85 .85 .85], 'LineStyle', 'none');hold on;
            uistack(p,'bottom')

            if i==1
                set(gca,'ytick',[0 1],'YTickLabels',[0 1], 'FontSize',7)

            else
                yticks([])
                yticklabels([])  
            end

            if j==axgrid(1) & i==1
                set(gca,'xtick',[0 .2],'XTickLabels',[0 200], 'FontSize',7)
                xtickangle(0)
            else
                xticks([])
                xtickangle(0)
                xticklabels([])
            end

         end


    end
    
end

