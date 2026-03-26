function good_chann = my_good_chann_TCI(indiv, sessN, dataFolder, binWidthMs, args)
% select good channels based on STRF, test retest on TCI & zeta test
arguments 
    indiv char
    sessN (1,1) double
    dataFolder char
    binWidthMs (1,1) double
    args.th (1,1) double = 0.1
    args.signalTime char = 'stim-aligned-v2'
    args.signalType char = 'quicksort/comavgrem'
    args.trialType char = 'passive'
    args.significanceMethod char = 'montecarlo'

end
trialType = args.trialType;
signalType = args.signalType;
signalTime = args.signalTime;
significanceMethod = args.significanceMethod;
th = args.th;

% load zeta
zetaParamString = 'dblUseMaxDur-0.75_intResampNum-1000_vecRestrictRange-0-0.5_boolDirectQuantile-1';

fname = strcat([ dataFolder filesep 'sub-' indiv filesep 'sess-'  sprintf('%03d',sessN) filesep 'puretones' filesep trialType filesep ... 
                                'preprocessing' filesep 'continuous' filesep 'quicksort/comavgrem' ]);


load(strcat([fname filesep zetaParamString filesep 'zeta.mat']));
good_chann_zeta = find([zeta.dblZetaP] <= 0.01);


%%
if strcmp(indiv, 'oscipek') | strcmp(indiv, 'wimereux') | (strcmp(indiv, 'murols') & sessN==35)
    STRF_bad_channels_mask_v2 = zeros(1,numel([zeta.dblZetaP]));
else
       
    strfParam = 'lags-0-50/method-ridge_nfolds-5_lagms-0-250_nlags-50';

    fname = strcat([ dataFolder filesep 'sub-' lower(indiv) filesep 'sess-' sprintf('%03d',sessN) filesep 'torcs' filesep trialType filesep ... 
                'STRF' filesep 'stim-aligned' filesep 'quicksort/comavgrem'  filesep ... 
                'binWidth-' num2str(5) 'ms' filesep strfParam]);

    load(strcat([fname filesep 'STRF_bad_channels_mask_v2.mat']));
end

%%
modelWinparam = 'model_fit_lossfn-sqerr_shape-1-2-3-4-5_intper_range-0.0078125-0.5_delay_range-0-0.25_nintper-100';
modelFitparam = 'boundstrength-0-0.25-0.5-1-2_nullsmps-100.mat';

taskType = 'contstream';

fileName = strcat([ dataFolder filesep 'sub-' indiv filesep 'sess-' sprintf('%03d',sessN) filesep taskType filesep trialType filesep ... 
                    'TCI' filesep signalTime filesep signalType filesep 'binWidth-' num2str(binWidthMs) 'ms' ...
                     filesep modelWinparam filesep modelFitparam]); % '.mat']);% 

%LOAD m
load(strcat([fileName]));

good_chann_tci = squeeze(median(M.same_context(:,1,:),1))'>th;
% -log10(.001) = 3 so that's what we threshold at
switch significanceMethod
    case 'guass_fit'
        fit_significance = M.logP_gaussfit'<= 0.001; % p=0.05
    case 'montecarlo'
        fit_significance = M.logP_counts'>=2 ;
    otherwise
        error('no matching significance testing method')
end

% good_chann_tci = any(squeeze(median(M.same_context,1))>th,1);


% good_chann = intersect(intersect(good_chann_zeta, good_chann_tci),find(~STRF_bad_channels_mask_v & 2));
% good_chann = intersect(good_chann_tci, find(~STRF_bad_channels_mask_v2));
good_chann = find(good_chann_tci & ~STRF_bad_channels_mask_v2 & fit_significance);
% good_chann = good_chann_zeta;



end