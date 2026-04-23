% Loads published data and reproduces plots.
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
figFolder = strcat([ 'figs' ]);

subnames={'A','B'};
sessions = {[31:40], ...
    [37,38,41,42]};

% for plotting 
markerstyles = {'^','square'}; 
cols = lines(7);
cols(2,:)=[];

rates = {'normal','compressed','stretched'};
categories = {'ferrets','speech'};
context={'random-random', 'random-natural'};



%% The effect of rate & category ->  y ∼1 + rate + category + (1 | unit) + (1 | session)

xi=1;
xii=1;

all_widths = zeros(1,500);
all_centers = zeros(1,500);
all_widths_wide = zeros(500,length(rates),length(categories),2,2);
all_distances = zeros(1,500);
all_subs = zeros(1,500);
all_subs_wide = zeros(1,500);
all_areas = cell(1,500);
all_unique_sess = cell(1,500);
all_unique_unit = cell(1,500);
all_recording = cell(1,500);
all_rates = cell(1,500);
all_cats = cell(1,500);

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

        all_subs_wide(xii:xii+units_n-1)=subid;
        % pickup model fits and distance
        for r=1:length(rates)
            for c=1:length(categories)
            
                all_widths(xi:xi+units_n-1)=M{r,c}.best_intper_sec;
                all_centers(xi:xi+units_n-1)=M{r,c}.best_delay_sec_median;
                all_distances(xi:xi+units_n-1)=distances;
                all_subs(xi:xi+units_n-1)=subid;
                all_areas(xi:xi+units_n-1)=areas;

                all_unique_sess(xi:xi+units_n-1) = repmat( {[sprintf('%03d',sessN) subjid ]}, 1, units_n);
                all_unique_unit(xi:xi+units_n-1) = arrayfun(@(x) [sprintf('%03d',sessN) subjid num2str(x) ], 1:units_n, 'UniformOutput', false);
                all_cats(xi:xi+units_n-1)= categories(c);
                all_rates(xi:xi+units_n-1)= rates(r);
                
                all_widths_wide(xii:xii+units_n-1,r,c,:,:)=M{r,c}.splits_best_intper_sec;
                
                
                xi=xi+units_n;
            end
            
        end
        all_recording(xii:xii+units_n-1)=repmat( {[sprintf('%03d',sessN) subjid ]}, 1, units_n);
     
        xii=xii+units_n;
    end
end


all_widths = all_widths(1:xi-1); all_widths=all_widths*1000;
all_centers = all_centers(1:xi-1); all_centers=all_centers*1000;
all_distances = all_distances(1:xi-1);
all_areas = all_areas(1:xi-1);
all_subs = all_subs(1:xi-1);
all_unique_sess = all_unique_sess(1:xi-1);
all_unique_unit = all_unique_unit(1:xi-1);
all_cats = all_cats(1:xi-1);
all_rates = all_rates(1:xi-1);

all_widths_wide = all_widths_wide(1:xii-1,:,:,:,:); all_widths_wide=all_widths_wide*1000;
all_subs_wide = all_subs_wide(1:xii-1); 
all_recording=all_recording(1:xii-1);
%% LME 
tbl = table(log2(all_widths(:)),all_rates(:),all_cats(:),all_unique_unit(:),all_unique_sess(:),'VariableNames',{'width','rate','category','unit','recording'});
lme = fitlme(tbl,'width ~ 1 + rate + category + (1 | unit) + (1 | recording)', 'CovariancePattern', {'diagonal','diagonal'})
fixedEffects = anova(lme,'DFMethod', 'satterthwaite')

fprintf(['rate: (F_{%d,%.2f} = %.2f, p = %.4f, \\beta_{%s} = %.3f octaves/mm, CI = [%.3f %.3f], \n' ...
    '\\beta_{%s} = %.3f octaves/mm, CI = [%.3f %.3f]; N = %d [observations for %d neurons])\n'],...,
    2,fixedEffects.DF2(2), fixedEffects.FStat(2), fixedEffects.pValue(2), lme.CoefficientNames{2}, ...
    lme.Coefficients.Estimate(2), lme.Coefficients.Lower(2), lme.Coefficients.Upper(2), lme.CoefficientNames{3}, ...
    lme.Coefficients.Estimate(3), lme.Coefficients.Lower(3), lme.Coefficients.Upper(3), ...
    lme.NumObservations,numel(lme.VariableInfo('unit',:).Range{1})) ;

%% bootrapped rank correlations, splits across stimuli
%%
% average widths across stimuli and rates for two data splits
all_widths_avgrate_splits =  squeeze(mean(all_widths_wide(:,:,:,[1,2],2),[3,2]));

statistic = @(A)rankcorr(A(:,1),A(:,2));
samples_overall = hierarch_bstrap(all_widths_avgrate_splits, all_recording, statistic, 10000);
p=mean(samples_overall<0)*2;
if p==0;p=1/numel(samples_overall);end

rho=rankcorr(all_widths_avgrate_splits(:,1),all_widths_avgrate_splits(:,2));

disp(strcat(['rank correlation, widths across splits, rho = ' num2str(rho) ', p=' num2str(p)]));

figure('Units','pixels','Position',[50 50 300 300]);
ax1 = axes(gcf,'Units','pixels','Position',[70 60 200 200 ]);  hold on
axis(gca, 'square'); box off;
hold on
xi = linspace(1, 500, 1000);

plot( xi, xi, '-','Linewidth',.5, 'Color',[.8 .8 .8])
scatter(all_widths_avgrate_splits(all_subs_wide==1,1),all_widths_avgrate_splits(all_subs_wide==1,2),10,'marker',markerstyles{1}, ...
            'MarkerEdgeColor', [.3 .3 .3],...
            'MarkerEdgeAlpha' ,1,...
            'MarkerFaceAlpha',.2,...
            'MarkerFaceColor', [.3 .3 .3]);
scatter(all_widths_avgrate_splits(all_subs_wide==2,1),all_widths_avgrate_splits(all_subs_wide==2,2),10,'marker',markerstyles{2}, ...
            'MarkerEdgeColor', [.3 .3 .3],...
            'MarkerEdgeAlpha' ,1,...
            'MarkerFaceAlpha',.2,...
            'MarkerFaceColor', [.3 .3 .3]);

set(gca,'YScale','log')
set(gca,'XScale','log')

ticks = [15 32 64 128 256];
xlim([12 220])
ylim([12 220])

set(gca, 'xtick', ticks, 'XTickLabels', ticks,'XMinorTick','off')
set(gca, 'Ytick', ticks, 'YTickLabels', ticks,'YMinorTick','off')



title({['avg window'],['rho = ' num2str(rho) ', p=' num2str(p)] ,[]});


xlabel(['Integraiton window (Split 1, ms)'], 'interpreter', 'none');
ylabel(['Integraiton window (Split 2, ms)'], 'interpreter', 'none');

set(findall(gcf,'-property','FontSize'),'FontSize',9)
exportgraphics(gcf, [figFolder filesep 'fig4e.png'] ,'Resolution',2000,'BackgroundColor','none')









%%
% collect the yoking metric for two data splits
yoking_metric_splits=[];

% average across category
% pickup widths for compressed and stretched,split across stimuli, average across categories
f=squeeze(mean(all_widths_wide(:,[2,3],:,[1,2],2),3));
yoking_metric_splits(:,1) = myYokingMetric(f(:,1,1),f(:,2,1),2);
yoking_metric_splits(:,2) = myYokingMetric(f(:,1,2),f(:,2,2),2);



samples_overall = hierarch_bstrap(yoking_metric_splits, all_recording, statistic, 10000);
p=mean(samples_overall<0)*2;
if p==0;p=1/numel(samples_overall);end

rho=rankcorr(yoking_metric_splits(:,1),yoking_metric_splits(:,2));



disp(strcat(['rank correlation, yoking metric across splits, rho = ' num2str(rho) ', p=' num2str(p)]));


figure('Units','pixels','Position',[50 50 300 300]);
ax1 = axes(gcf,'Units','pixels','Position',[70 60 200 200 ]);  hold on
axis(gca, 'square'); box off;
hold on


plot([-2 2],[-2 2], '-','Linewidth',.5, 'Color',[.8 .8 .8])
 xline(0,'color', cols(4,:),'Linewidth',1.5,'Layer','bottom')

yline(0,'color', cols(4,:),'Linewidth',1.5,'Layer','bottom')

xline(1,'color', cols(3,:),'Linewidth',1.5,'Layer','bottom')
yline(1,'color', cols(3,:),'Linewidth',1.5,'Layer','bottom')


scatter(yoking_metric_splits(all_subs_wide==1,1),yoking_metric_splits(all_subs_wide==1,2),10,'marker',markerstyles{1}, ...
            'MarkerEdgeColor', [.3 .3 .3],...
            'MarkerEdgeAlpha' ,1,...
            'MarkerFaceAlpha',.2,...
            'MarkerFaceColor', [.3 .3 .3]);
   
scatter(yoking_metric_splits(all_subs_wide==2,1),yoking_metric_splits(all_subs_wide==2,2),10,'marker',markerstyles{2}, ...
            'MarkerEdgeColor', [.3 .3 .3],...
            'MarkerEdgeAlpha' ,1,...
            'MarkerFaceAlpha',.2,...
            'MarkerFaceColor', [.3 .3 .3]);
  


title({['rate-yoking index'],['rho = ' num2str(rho) ', p=' num2str(p)],[] });

xlim([-1.9 1.9])
ylim([-1.9 1.9])
xticks([-1 0 1])
yticks([-1 0 1])

xlabel(['Rate-yoking index (Split 1)'], 'interpreter', 'none');
ylabel(['Rate-yoking index (Split 2)'], 'interpreter', 'none');
exportgraphics(gcf, [figFolder filesep 'fig4d.png'] ,'Resolution',2000,'BackgroundColor','none')

%% scatter 
   
figure('Units','pixels','Position',[50 50 300 300]);
ax1 = axes(gcf,'Units','pixels','Position',[70 70 200 200 ]);  hold on
axis(gca, 'square'); box off;

hold on
myFormatScatter(gca,sqrt(2)/(1/sqrt(2)))
for subid= 1:length(subnames)

 
    markerstyle = markerstyles{subid};
    sh=[];
    for split=1:2
        stretched = (all_widths_wide(all_subs_wide==subid,3,1,split,2) + all_widths_wide(all_subs_wide==subid,3,2,split,2))/2 ;
     
        compressed = (all_widths_wide(all_subs_wide==subid,2,1,split,2) + all_widths_wide(all_subs_wide==subid,2,2,split,2))/2 ;
 
        
        sh{split} = scatter(stretched,compressed,10,'marker',markerstyle, ...
                    'MarkerEdgeColor', [.3 .3 .3],...
                    'MarkerEdgeAlpha' ,.7,...
                    'MarkerFaceAlpha',.2,...
                    'MarkerFaceColor', [.3 .3 .3]);
    end

    plot([sh{1}.XData' sh{2}.XData' ]',[sh{1}.YData'  sh{2}.YData']','Color',[.8 .8 .8 .5])

end

ylabel([rates{2} ' (ms)'], 'interpreter', 'none');
xlabel([rates{3} ' (ms)'], 'interpreter', 'none');
set(findall(gcf,'-property','FontSize'),'FontSize',14)
exportgraphics(gcf, [figFolder filesep 'fig4c.png'] ,'Resolution',2000,'BackgroundColor','none')


%% scatter catgories

cols=lines(3);cols=cols([1,3],:)
figure('Units','pixels','Position',[50 50 300 300]);
ax1 = axes(gcf,'Units','pixels','Position',[70 70 200 200 ]);  hold on
axis(gca, 'square'); box off;

hold on
myFormatScatter(gca,sqrt(2)/(1/sqrt(2)))
hold on
for subid= 1:length(subnames)

 
    markerstyle = markerstyles{subid};
    sh=[];
    for split=1:2

      

        % both categories

        stretched = squeeze(all_widths_wide(all_subs_wide==subid,3,:,split,2));
     
        compressed = squeeze(all_widths_wide(all_subs_wide==subid,2,:,split,2)) ;

  
         
        sh{1,split} = scatter(stretched(:,1),compressed(:,1),7,'marker',markerstyle, ...
                    'MarkerEdgeColor',cols(1,:),...
                    'MarkerEdgeAlpha' ,.5,...
                    'MarkerFaceAlpha',.1,...
                    'MarkerFaceColor',cols(1,:),...
                    'jitter','on','JitterAmount',.3,'DisplayName','ferret');
        sh{2,split} =scatter(stretched(:,2),compressed(:,2),7,'marker',markerstyle, ...
                    'MarkerEdgeColor',cols(2,:),...
                    'MarkerEdgeAlpha' ,.5,...
                    'MarkerFaceAlpha',.1,...
                    'MarkerFaceColor',cols(2,:),...
                    'jitter','on','JitterAmount',.3,'DisplayName','speech');
    end
    plot([[sh{1,1}.XData]' [sh{1,2}.XData]' ]',[[sh{1,1}.YData]'  [sh{1,2}.YData]']','Color', [cols(1,:) .3])
    plot([[sh{2,1}.XData]' [sh{2,2}.XData]' ]',[[sh{2,1}.YData]'  [sh{2,2}.YData]']','Color', [cols(2,:) .3])
end
    
xlim([12 250])
ylim([12 250])


ylabel([rates{2} ' (ms)'], 'interpreter', 'none');
xlabel([rates{3} ' (ms)'], 'interpreter', 'none');

set(findall(gcf,'-property','FontSize'),'FontSize',14)
exportgraphics(gcf, [figFolder filesep 'fig4c_catsup.png'] ,'Resolution',2000,'BackgroundColor','none')


%% identity scatter per category
    
figure('Units','pixels','Position',[50 50 400 400]);
ax1 = axes(gcf,'Units','pixels','Position',[70 60 200 200 ]);  hold on
axis(gca, 'square'); box off;
hold on
xi = linspace(1, 500, 1000);

plot( xi, xi, '-','Linewidth',.5, 'Color',[.8 .8 .8])
for c=1:2

    % 2categore across, axis is splits
   

    x = squeeze(mean(all_widths_wide(:,:,c,:,2),2));
 
    speech = squeeze(all_widths_wide(all_subs_wide==subid,2,:,split,2)) ;

    scatter(x(all_subs_wide==1,1),x(all_subs_wide==1,2),10,'marker',markerstyles{1}, ...
                'MarkerEdgeColor', cols(c,:),...
                'MarkerEdgeAlpha' ,.5,...
                'MarkerFaceAlpha',.2,...
                'MarkerFaceColor', cols(c,:));
    scatter(x(all_subs_wide==2,1),x(all_subs_wide==2,2),10,'marker',markerstyles{2}, ...
                'MarkerEdgeColor', cols(c,:),...
                'MarkerEdgeAlpha' ,.5,...
                'MarkerFaceAlpha',.2,...
                'MarkerFaceColor', cols(c,:));
end

set(gca,'YScale','log')
set(gca,'XScale','log')

ticks = [15 32 64 128 256];
xlim([12 220])
ylim([12 220])

set(gca, 'xtick', ticks, 'XTickLabels', ticks,'XMinorTick','off')
set(gca, 'Ytick', ticks, 'YTickLabels', ticks,'YMinorTick','off')


xlabel(['Integraiton window (Split 1, ms)'], 'interpreter', 'none');
ylabel(['Integraiton window (Split 2, ms)'], 'interpreter', 'none');

set(findall(gcf,'-property','FontSize'),'FontSize',9)
exportgraphics(gcf, [figFolder filesep 'fig4e_catsup.png'] ,'Resolution',2000,'BackgroundColor','none')

%% yoking metric per scatter 
cols = lines(7);
cols(2,:)=[];
% collect the yoking metric for two data splits
yoking_metric_splits=[];

% per category
% pickup widths for compressed and stretched,split across stimuli, average across categories
f=squeeze(all_widths_wide(:,[2,3],:,[1,2],2));
%ferret
yoking_metric_splits_cat(1,:,1) = myYokingMetric(f(:,1,1,1),f(:,2,1,1),2);
yoking_metric_splits_cat(1,:,2) = myYokingMetric(f(:,1,1,2),f(:,2,1,2),2);
%speech
yoking_metric_splits_cat(2,:,1) = myYokingMetric(f(:,1,2,1),f(:,2,2,1),2);
yoking_metric_splits_cat(2,:,2) = myYokingMetric(f(:,1,2,2),f(:,2,2,2),2);

figure('Units','pixels','Position',[50 50 300 300]);
ax1 = axes(gcf,'Units','pixels','Position',[70 60 200 200 ]);  hold on
axis(gca, 'square'); box off;
hold on


plot([-2 2],[-2 2], '-','Linewidth',.5, 'Color',[.8 .8 .8])
 xline(0,'color', cols(4,:),'Linewidth',1.5,'Layer','bottom')

yline(0,'color', cols(4,:),'Linewidth',1.5,'Layer','bottom')

xline(1,'color', cols(3,:),'Linewidth',1.5,'Layer','bottom')
yline(1,'color', cols(3,:),'Linewidth',1.5,'Layer','bottom')


scatter(yoking_metric_splits_cat(1,all_subs_wide==1,1),yoking_metric_splits_cat(1,all_subs_wide==1,2),10,'marker',markerstyles{1}, ...
            'MarkerEdgeColor', cols(1,:),...
            'MarkerEdgeAlpha' ,1,...
            'MarkerFaceAlpha',.2,...
            'MarkerFaceColor', cols(1,:));
   
scatter(yoking_metric_splits_cat(1,all_subs_wide==2,1),yoking_metric_splits_cat(1,all_subs_wide==2,2),10,'marker',markerstyles{2}, ...
            'MarkerEdgeColor', cols(1,:),...
            'MarkerEdgeAlpha' ,1,...
            'MarkerFaceAlpha',.2,...
            'MarkerFaceColor',cols(1,:));
  

scatter(yoking_metric_splits_cat(2,all_subs_wide==1,1),yoking_metric_splits_cat(2,all_subs_wide==1,2),10,'marker',markerstyles{1}, ...
            'MarkerEdgeColor', cols(2,:),...
            'MarkerEdgeAlpha' ,1,...
            'MarkerFaceAlpha',.2,...
            'MarkerFaceColor', cols(2,:));
   
scatter(yoking_metric_splits_cat(2,all_subs_wide==2,1),yoking_metric_splits_cat(2,all_subs_wide==2,2),10,'marker',markerstyles{2}, ...
            'MarkerEdgeColor', cols(2,:),...
            'MarkerEdgeAlpha' ,1,...
            'MarkerFaceAlpha',.2,...
            'MarkerFaceColor',cols(2,:));


xlim([-1.9 1.9])
ylim([-1.9 1.9])
xticks([-1 0 1])
yticks([-1 0 1])

xlabel(['Rate-yoking index (Split 1)'], 'interpreter', 'none');
ylabel(['Rate-yoking index (Split 2)'], 'interpreter', 'none');
exportgraphics(gcf, [figFolder filesep 'fig4d_catsup.png'] ,'Resolution',2000,'BackgroundColor','none')





















function metric = myYokingMetric(x,y,factor)
     % x ~ compressed
     % y ~ stretched

    metric = (log2(x)-log2(y))/log2(factor);

    

end


function myFormatScatter(ax,factor)
    cols = lines(5);
    

    xi = linspace(1, 500, 1000);
    if factor==1
        plot(ax,xi, xi, '-','Linewidth',1.5, 'Color', cols(5,:)); hold on;
    else
        
        plot(ax,xi, xi, '-','Linewidth',1.5, 'Color', cols(5,:)); hold on;
    end
    plot(ax,  xi*factor,xi, '-','Linewidth',1.5, 'Color',cols(4,:))
    set(ax,'YScale','log')
    set(ax,'XScale','log')

    ticks = [16 32 64 128 256];
    xlim([12 220])
    ylim([12 220])

    set(ax, 'xtick', ticks, 'XTickLabels', ticks,'XMinorTick','off')
    set(ax, 'Ytick', ticks, 'YTickLabels', ticks,'YMinorTick','off')


end
