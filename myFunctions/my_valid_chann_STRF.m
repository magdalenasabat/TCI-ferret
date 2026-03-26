function good_chann = my_valid_chann_STRF(indiv, sessN, dataFolder,  args)
% select valid channels based on STRF
% some channles chave sound artifacts, we catch it based on instantenaus
% response in STRF
arguments 
    indiv char
    sessN (1,1) double
    dataFolder char
    args.signalType char = 'quicksort/comavgrem'
    args.trialType char = 'passive'

end
trialType = args.trialType;
signalType = args.signalType;

%%
if strcmp(indiv, 'oscipek') | strcmp(indiv, 'wimereux') | (strcmp(indiv, 'murols') & sessN==35)
    fprintf('No STRF - all valid chans')
    STRF_bad_channels_mask_v2 = false;
else
       
    strfParam = 'lags-0-50/method-ridge_nfolds-5_lagms-0-250_nlags-50';

    fname = strcat([ dataFolder filesep 'sub-' lower(indiv) filesep 'sess-' sprintf('%03d',sessN) filesep 'torcs' filesep trialType filesep ... 
                'STRF' filesep 'stim-aligned' filesep signalType  filesep ... 
                'binWidth-' num2str(5) 'ms' filesep strfParam]);

    load(strcat([fname filesep 'STRF_bad_channels_mask_v2.mat']));
end

%%

good_chann =  ~STRF_bad_channels_mask_v2 ;


end