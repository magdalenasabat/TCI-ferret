function good_chann = my_good_chann_TCI_slow(indiv, sessN, dataFolder, binWidthMs, args)
% select good channels based on STRF, test retest on TCI & zeta test
arguments 
    indiv char
    sessN (1,1) double
    dataFolder char
    binWidthMs (1,1) double
    args.th (1,1) double = 0.1
    args.signalTime char = 'stim-aligned-v2'
    args.signalType char = 'quicksort/comavgrem'
    args.taskType char = 'tci-slow';
    args.trialType char = 'passive'

end
trialType = args.trialType;
signalType = args.signalType;
signalTime = args.signalTime;
taskType = args.taskType;
th = args.th;

%%


fileName = strcat([ dataFolder filesep 'sub-' indiv filesep 'sess-' sprintf('%03d',sessN) filesep taskType filesep trialType filesep ... 
                    'TCI' filesep signalTime filesep signalType filesep 'binWidth-' num2str(binWidthMs) 'ms' ...
                     filesep  filesep 'sub-' indiv '_ses-' sprintf('%03d',sessN) '_task-' taskType ' _trial-' trialType '_Lstruct.mat' ]);

%LOAD m
load(strcat([fileName]));

good_chann_tci = squeeze(median(L.pooled.same_context(:,1,:),1))'>th;
% good_chann_tci = any(squeeze(median(M.same_context,1))>th,1);

%%
if strcmp(indiv, 'oscipek') | strcmp(indiv, 'wimereux') | ((strcmp(indiv, 'murols') & sessN==35))
    STRF_bad_channels_mask_v2 = zeros(1,length(good_chann_tci));
else
       
    strfParam = 'lags-0-50/method-ridge_nfolds-5_lagms-0-250_nlags-50';

    fname = strcat([ dataFolder filesep 'sub-' lower(indiv) filesep 'sess-' sprintf('%03d',sessN) filesep 'torcs' filesep trialType filesep ... 
                'STRF' filesep 'stim-aligned' filesep 'quicksort/comavgrem'  filesep ... 
                'binWidth-' num2str(5) 'ms' filesep strfParam]);

    load(strcat([fname filesep 'STRF_bad_channels_mask_v2.mat']));
end




% good_chann = intersect(intersect(good_chann_zeta, good_chann_tci),find(~STRF_bad_channels_mask_v2));
% good_chann = intersect(good_chann_tci, find(~STRF_bad_channels_mask_v2));
good_chann = find(good_chann_tci & ~STRF_bad_channels_mask_v2);
% good_chann = good_chann_zeta;



end