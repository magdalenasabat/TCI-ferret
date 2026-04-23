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
context={'random-random', 'random-natural'};

%% The effect of rate & category ->  y ∼1 + rate + category + (1 | unit) + (1 | session)

xi=1;
all_widths = zeros(1,500);
all_widths_splits = zeros(2,500);
all_centers = zeros(1,500);
all_distances = zeros(1,500);
all_subs = zeros(1,500);
all_areas = cell(1,500);
all_unique_sess = cell(1,500);
all_unique_unit = cell(1,500);
all_rates = cell(1,500);
all_context = cell(1,500);

for subid= 1:length(subnames) 
    
    subjid = subnames{subid};
    fprintf(strcat([' SUBJECT ' subjid '\n']))
    
    
    for sessN=sessions{subid}
        
%         fprintf(strcat([' sess #' sprintf('%03d',sessN) '\n']))

        %load files
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_anatomical.mat'  ]));
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_model_fits.mat'  ]));
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_cross_context_correlation.mat'  ]));
 
        units_n = length(M_context{1}.best_intper_sec);
        
        % pickup model fits and distance
        for r=1:length(rates)
            for c=1:length(context)
            
                all_widths(xi:xi+units_n-1)=M_context{r,c}.best_intper_sec;
                if units_n>1
                    all_widths_splits(:,xi:xi+units_n-1)=M_context{r,c}.splits_best_intper_sec(:,:,2)';
                else
                    all_widths_splits(:,xi:xi+units_n-1)=M_context{r,c}.splits_best_intper_sec(:,2)';
                end
                all_centers(xi:xi+units_n-1)=M_context{r,c}.best_delay_sec_median;
                all_distances(xi:xi+units_n-1)=distances;
                all_subs(xi:xi+units_n-1)=subid;
                all_areas(xi:xi+units_n-1)=areas;

                all_unique_sess(xi:xi+units_n-1) = repmat( {[sprintf('%03d',sessN) subjid ]}, 1, units_n);
                all_unique_unit(xi:xi+units_n-1) = arrayfun(@(x) [sprintf('%03d',sessN) subjid num2str(x) ], 1:units_n, 'UniformOutput', false);
                all_context(xi:xi+units_n-1)= context(c);
                all_rates(xi:xi+units_n-1)= rates(r);
                xi=xi+units_n;
            end
        end
        
    end
end


all_widths = all_widths(1:xi-1); all_widths=all_widths*1000;
all_widths_splits = all_widths_splits(:,1:xi-1); all_widths_splits=all_widths_splits*1000;
all_centers = all_centers(1:xi-1); all_centers=all_centers*1000;
all_distances = all_distances(1:xi-1);
all_areas = all_areas(1:xi-1);
all_subs = all_subs(1:xi-1);
all_unique_sess = all_unique_sess(1:xi-1);
all_unique_unit = all_unique_unit(1:xi-1);
all_context = all_context(1:xi-1);
all_rates = all_rates(1:xi-1);

%% LME 
tbl = table(log2(all_widths(:)),all_rates(:),all_context(:),all_unique_unit(:),all_unique_sess(:),'VariableNames',{'width','rate','context','unit','recording'});
lme = fitlme(tbl,'width ~ 1 + rate + context + (1 | unit) + (1 | recording)', 'CovariancePattern', {'diagonal','diagonal'})
fixedEffects = anova(lme,'DFMethod', 'satterthwaite')

fprintf(['rate: (F_{%d,%.2f} = %.2f, p = %.4f, \\beta_{%s} = %.3f octaves/mm, CI = [%.3f %.3f], \n' ...
    '\\beta_{%s} = %.3f octaves/mm, CI = [%.3f %.3f]; N = %d [observations for %d neurons])\n', ...
    'context: (F_{%d,%.2f} = %.2f, p = %.4f, \\beta_{%s} = %.3f octaves, CI = [%.3f %.3f]; \n' ...
    ' N = %d [observations for %d neurons])\n'],...,
    2,fixedEffects.DF2(2), fixedEffects.FStat(2), fixedEffects.pValue(2), lme.CoefficientNames{2}, ...
    lme.Coefficients.Estimate(2), lme.Coefficients.Lower(2), lme.Coefficients.Upper(2), lme.CoefficientNames{3}, ...
    lme.Coefficients.Estimate(3), lme.Coefficients.Lower(3), lme.Coefficients.Upper(3), ...
    lme.NumObservations,numel(lme.VariableInfo('unit',:).Range{1}), ...
    1,fixedEffects.DF2(3), fixedEffects.FStat(3), fixedEffects.pValue(3), lme.CoefficientNames{4}, ...
    lme.Coefficients.Estimate(4), lme.Coefficients.Lower(4), lme.Coefficients.Upper(4), ...
    lme.NumObservations,numel(lme.VariableInfo('unit',:).Range{1}) );
%%
addpath(genpath(strcat(['code' filesep 'toolbox-ScientificColourMaps8'])));

load('roma.mat');
hist_col = roma([70,10,230],:);
bins = 10.^(1:.11:2.4);


for c=unique(all_context)
    f = figure('Units', 'centimeters','Position', [50,50,2.5,2.15]);
    hold on;
    ri=1;
    for r=unique(all_rates)


        h=histogram(all_widths(strcmp(all_context,c) & strcmp(all_rates,r) ),bins,...
            'FaceColor',hist_col(ri,:),'EdgeColor',hist_col(ri,:),'FaceAlpha',.7); %',.4,[.5 .5 .5],lag_t(2:end-1),[], 'sem');
        
        set(gca,'XScale','log')
        ticks = [16  32 64 128 ];
        xlim([min(ticks)-6 max(ticks)+64])
        xtickangle(0)
        if ri ==3
            set(gca, 'xtick', ticks, 'XTickLabels', ticks,'XMinorTick','off')
        else
            set(gca, 'xtick', ticks, 'XTickLabels', [],'XMinorTick','off')
        end
        set(gca, 'ytick', [0:20:100])
        ri=ri+1;
    end
    set(findall(f,'-property','FontSize'),'FontSize',7)
    exportgraphics(gcf, [dataFolder filesep  c{:} '_hist.png'] ,'Resolution',2000,'BackgroundColor','none')
    
end


%%

A = all_widths_splits(:,strcmp(all_context,'random-random'))-all_widths_splits(:,strcmp(all_context,'random-natural'));

statistic = @(A)rankcorr(A(:,1),A(:,2));
samples_overall = hierarch_bstrap(A', all_unique_sess(strcmp(all_context,'random-random'))', statistic, 10000);
p=mean(samples_overall<0)*2;
if p==0;p=1/numel(samples_overall);end

rho=rankcorr(A(1,:)',A(2,:)');

disp(strcat(['rank correlation, context across splits, rho = ' num2str(rho) ', p=' num2str(p)]));