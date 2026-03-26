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

dataFolder = strcat([ 'data' filesep 'experiment1']);

subnames={'A','B','C'};
sessions = {[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,26,27,29,30,31], ...
    [1,2,3,4,5,7,8,9,10,12,13,14,15,22,23,25,26,27,28,30,31,32,34,35,38],...
    [2,3,4,5,6,7,9,10,11,12,13,14,15]};

% for plotting 
markerstyles = {'^','square', 'o'}; 

%% get data

xi=1;
all_widths = zeros(1,500);
all_centers = zeros(1,500);
all_distances = zeros(1,500);
all_subs = zeros(1,500);
all_areas = cell(1,500);
all_unique_sess = cell(1,500);
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
        all_widths(xi:xi+units_n-1)=M.best_intper_sec;
        all_centers(xi:xi+units_n-1)=M.best_delay_sec_median;
        all_distances(xi:xi+units_n-1)=distances;
        all_subs(xi:xi+units_n-1)=subid;
        all_areas(xi:xi+units_n-1)=areas;
        
        all_CCC(:,:,xi:xi+units_n-1)=L.diff_context;
        all_NC(:,:,xi:xi+units_n-1)=L.same_context;
        
        
        all_unique_sess(xi:xi+units_n-1) = repmat( {[sprintf('%03d',sessN) subjid ]}, 1, units_n);
        xi=xi+units_n;
        
    end
end


all_widths = all_widths(1:xi-1); all_widths=all_widths*1000;
all_centers = all_centers(1:xi-1); all_centers=all_centers*1000;
all_distances = all_distances(1:xi-1);
all_areas = all_areas(1:xi-1);
all_subs = all_subs(1:xi-1);
all_CCC = all_CCC(:,:,1:xi-1);
all_NC = all_NC(:,:,1:xi-1);
all_unique_sess = all_unique_sess(1:xi-1);
lag_t = L.lag_t;
unique_segs = L.unique_segs;

%% LME region and distance effects

fprintf(['Median width in the primary auditory cortex: %.1f ms, %.2f ms, N = %.0f \n',...
    'Median width in the non-primary auditory cortex: %.1f ms, %.2f ms,N = %.0f\n'], ...
    median(all_widths(strcmp(all_areas,'A1'))), std(all_widths(strcmp(all_areas,'A1'))), numel(all_widths(strcmp(all_areas,'A1'))),...
    median(all_widths(strcmp(all_areas,'PEG'))), std(all_widths(strcmp(all_areas,'PEG'))), numel(all_widths(strcmp(all_areas,'PEG'))) ...
    )


tbl = table(log2(all_widths(:)),log2(all_centers(:)),all_distances(:),all_areas(:),all_unique_sess(:),'VariableNames',{'width','center','distance','region','recording'});


% EFFECT OF DISTANCE TO PAC
lme = fitlme(tbl,'width ~ 1 + distance  + (1 + distance|recording)','CovariancePattern','diagonal');
fixedEffects = anova(lme,'DFMethod', 'satterthwaite');
fprintf(['(F_{%d,%.2f} = %.2f, p = %s, \\beta_{%s} = %.2f octaves/mm, CI = [%.3f %.3f]; \n' ...
    ' N = %d)\n'], ...
    1,fixedEffects.DF2(2), fixedEffects.FStat(2), fixedEffects.pValue(2), lme.CoefficientNames{2}, ...
    lme.Coefficients.Estimate(2), lme.Coefficients.Lower(2), lme.Coefficients.Upper(2), ...
    lme.NumObservations);
% EFFECT OF REGION
lme = fitlme(tbl,'width ~ 1 + region + (1 +region | recording)','CovariancePattern','diagonal');
fixedEffects = anova(lme,'DFMethod', 'satterthwaite');
fprintf(['(F_{%d,%.2f} = %.2f, p = %s, \\beta_{%s} = %.2f octaves/mm, CI = [%.3f %.3f]; \n' ...
    ' N = %d)\n'], ...
    1,fixedEffects.DF2(2), fixedEffects.FStat(2), fixedEffects.pValue(2), lme.CoefficientNames{2}, ...
    lme.Coefficients.Estimate(2), lme.Coefficients.Lower(2), lme.Coefficients.Upper(2), ...
    lme.NumObservations);

% fprintf(['Median center in the primary auditory cortex: %.2f ms, %.2f ms, N = %.0f \n',...
%     'Median center in the non-primary auditory cortex: %.2f ms, %.2f ms,'], ...
%     median(all_centers(strcmp(all_areas,'A1'))), std(all_centers(strcmp(all_areas,'A1'))), numel(all_centers(strcmp(all_areas,'A1'))),...
%     median(all_centers(strcmp(all_areas,'PEG'))), std(all_centers(strcmp(all_areas,'PEG'))), numel(all_centers(strcmp(all_areas,'PEG'))) ...
%     )
% lme = fitlme(tbl,'center ~ 1 + distance + (1+ distance|recording)','CovariancePattern','diagonal')
% lme = fitlme(tbl,'center ~ 1 + region + (1 +region | recording)','CovariancePattern','diagonal')
% 
