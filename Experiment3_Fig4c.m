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

dataFolder = strcat([ 'data' filesep 'experiment3']);

subnames={'A','B'};
sessions = {[31:40], ...
    [37,38,41,42]};

% for plotting 
markerstyles = {'^','square'}; 

rates = {'normal','compressed','stretched'};
categories = {'ferrets','speech'};
context={'random-random', 'random-natural'};

%% 

xi=1;
all_widths = zeros(1,500,length(rates),length(categories));
all_centers = zeros(1,500,length(rates),length(categories));
all_subs = zeros(1,500);

for subid= 1:length(subnames) 
    
    subjid = subnames{subid};
    fprintf(strcat([' SUBJECT ' subjid '\n']))
    
    
    for sessN=sessions{subid}
        
%         fprintf(strcat([' sess #' sprintf('%03d',sessN) '\n']))

        %load files
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_anatomical.mat'  ]));
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_model_fits.mat'  ]));
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_cross_context_correlation.mat'  ]));
 
        units_n = length(M{1}.best_intper_sec);
        
        % pickup model fits and distance
        for r=1:length(rates)
            for c=1:length(categories)
                n_seg= sum(L_context{r,c}.n_total_segs~=0);
                all_widths(xi:xi+units_n-1,r,c)=M{r,c}.best_intper_sec;
                all_centers(xi:xi+units_n-1,r,c)=M{r,c}.best_delay_sec_median;
                all_subs(xi:xi+units_n-1)=subid;

            end
        end
        xi=xi+units_n;
    end
end


all_widths = all_widths(1:xi-1,:,:); all_widths=all_widths*1000;
all_centers = all_centers(1:xi-1,:,:); all_centers=all_centers*1000;
all_subs = all_subs(1:xi-1);
%% Fig 4c


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

