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

%% The effect of rate & category ->  y ∼1 + rate + category + (1 | unit) + (1 | session)

xi=1;
xii=1;

all_widths = zeros(1,500);
all_centers = zeros(1,500);
all_widths_wide = zeros(500,length(rates),length(categories),2,2);
all_distances = zeros(1,500);
all_subs = zeros(1,500);
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

% average widths across stimuli and rates for two data splits
all_widths_avgrate_splits =  squeeze(mean(all_widths_wide(:,:,:,[1,2],2),[3,2]));

statistic = @(A)rankcorr(A(:,1),A(:,2));
samples_overall = hierarch_bstrap(all_widths_avgrate_splits, all_recording, statistic, 10000);
p=mean(samples_overall<0)*2;
if p==0;p=1/numel(samples_overall);end

rho=rankcorr(all_widths_avgrate_splits(:,1),all_widths_avgrate_splits(:,2));

disp(strcat(['rank correlation, widths across splits, rho = ' num2str(rho) ', p=' num2str(p)]));

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



function metric = myYokingMetric(x,y,factor)
     % x ~ compressed
     % y ~ stretched

    metric = (log2(x)-log2(y))/log2(factor);

    

end