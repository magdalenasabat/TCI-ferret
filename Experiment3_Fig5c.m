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

dataFolder = strcat([ 'data' filesep 'experiment3' ]);

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
all_CCC = NaN(696,6,500,length(rates),length(context));
all_NC = NaN(696,6,500,length(rates),length(context));

for subid= 1:length(subnames) 
    
    subjid = subnames{subid};
    fprintf(strcat([' SUBJECT ' subjid '\n']))
    
    
    for sessN=sessions{subid}
        
%         fprintf(strcat([' sess #' sprintf('%03d',sessN) '\n']))

        %load files
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_anatomical.mat'  ]));
       
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_cross_context_correlation.mat'  ]));
 
        units_n = length(L_context{1,1}.chnames);
        
        % pickup model fits and distance
        for r=1:length(rates)
            for c=1:length(context)
                n_seg= sum(L_context{r,c}.n_total_segs~=0);
               
                
                all_CCC(:,1:n_seg,xi:xi+units_n-1,r,c)=L_context{r,c}.diff_context(:,1:n_seg,:);
                all_NC(:,1:n_seg,xi:xi+units_n-1,r,c)=L_context{r,c}.same_context(:,1:n_seg,:);
        
             
            end
        end
        xi=xi+units_n;
    end
end


all_CCC = all_CCC(:,:,1:xi-1,:,:);
all_NC = all_NC(:,:,1:xi-1,:,:);

%% Fig 4a

addpath(genpath(strcat(['code' filesep 'toolbox-ScientificColourMaps8'])));
load('roma.mat');
hist_col = roma([70,10,230],:);
valid_durs = [1:4];
axgrid = [2,numel(valid_durs)];  % [#rows, #cols]
titles = {'closest', 'intermediate', 'farthest'}; 

figTCI = figure('Units', 'centimeters','Position', [50,50,7.9,5]);
% tclMain = tiledlayout(axgrid(1),1); 
% tcl = gobjects(1,axgrid(1));
% ax = gobjects(axgrid); 

titles={'Primary Auditory Cortex - MEG','Non-primary Auditory Cortex - PEG'};
idexes = [1:1:axgrid(2)];

axs = tiledlayout(axgrid(1),axgrid(2))

unique_segs=L_context{1, 1}.unique_segs;
lag_t=L_context{1, 2}.lag_t;


gr = unique(rates);
for c=1:axgrid(1)
    k=1;
    for i= 1:axgrid(2)
    %     if i==axgrid(2); continue; end
        idx = valid_durs(i);

        ax(c,k,i) = nexttile();
        hold( ax(c,k,i) , 'on')
        for j=1:length(gr)

            s = squeeze(all_NC(:,idx,:,j,c));
            d = squeeze(all_CCC(:,idx,:,j,c));

            plot(lag_t,nanmedian(d./s,2),'LineWidth',0.5,'color',hist_col(j,:))
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
    title([L_context{1,c}.context ' 5 ms'])

end

set(findall(figTCI,'-property','FontSize'),'FontSize',7)