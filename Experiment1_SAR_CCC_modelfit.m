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



addpath(genpath(strcat(['code' filesep 'toolbox-sam-tci'])));
addpath(genpath(strcat(['code' filesep 'myFunctions'])));

dataFolder = strcat([ 'data' filesep 'experiment1']);

subnames={'A','B','C'};


%% Example for the 1st session of animal A

subid = 1;
sessN = 1;
  
subjid = subnames{subid};
fprintf(strcat([' SUBJECT ' subjid '\n']))


load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_electrode_timecourses.mat'  ]));
load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_anatomical.mat'  ]));
load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_cross_context_correlation.mat'  ]));
load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_model_fits.mat'  ]));

%%
% % calculate correlations
% full output includes random-random comparison (L.samedur)
% random-natural comparison (L.diffdur) 
% and a pooled one (L.alldur) (all analysis exept Fig 5 where we compare the two context comparison types)
L_new = myCCC(c1,S,t,chnames,subjid,sessN,diffdur=c2 );

% compare to the published data
figure;
nexttile; plot(L_new.alldur.same_context(:,1,1)); hold on; plot(L_new.alldur.diff_context(:,1,1)); title('reproduced')
nexttile; plot(L.same_context(:,1,1)); hold on; plot(L.diff_context(:,1,1)); title('published')


%% fit the model
% in the published paper we tested a fairly wide parameter space,   
M_new = modelfit_cross_context_corr(L_new.alldur, 'overwrite', true, 'plot_figure', false, ...
                        'shape', [1,2,3,4,5], 'intper_range', [1/128, 0.50], 'delay_range', [0, 0.25], ...
                        'lossfn', 'sqerr', 'nintper', 100, 'boundstrength', [0, 0.25, 0.5, 1, 2]);

%%%%%%%%%%%%%%%%%
% % uncomment to use
% % alternatively try a smaller parameter set and compare model fits qualitatively 
%
% M_new = modelfit_cross_context_corr(L_new.alldur, 'overwrite', true, 'plot_figure', false, ...
%                         'shape', 1, 'intper_range', [1/128, 0.50], 'delay_range', [0, 0.25], ...
%                         'lossfn', 'sqerr');

% compare the fit to our data
figure;
nexttile; plot(L_new.alldur.same_context(:,1,1)); hold on; plot(L_new.alldur.diff_context(:,1,1)); plot(M_new.diff_context_bestpred(:,1,1));  
title({'reproduced', ['model width ' num2str(M_new.best_intper_sec(1))]})
nexttile; plot(L.same_context(:,1,1)); hold on; plot(L.diff_context(:,1,1));  plot(M.diff_context_bestpred(:,1,1)); 
title({'published', ['model width ' num2str(M.best_intper_sec(1))]})
legend('noise-ceiling','CCC','model prediction')






