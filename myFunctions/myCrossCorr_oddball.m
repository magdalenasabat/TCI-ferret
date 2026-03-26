function L = myCrossCorr_latest(dat, source_folder, plot_figure, smoothing )
signalType=dat.signalType;

if strcmp(signalType, 'spikeSorted')
    dat.processing_params = 'kilsort';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% intitialise

if ~exist('plot_figure','var')
  plot_figure = 'False';
end

figureFolder = strcat([ source_folder filesep 'analysis' filesep 'figures']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cleanup


fprintf('calculating correlations  \n')

% think about reorganising so that D A @ can be indexed
standard_same_context = NaN(size(dat.timebins,2),2, size(dat.binned,2));
standard_diff_context = NaN(size(dat.timebins,2), 2, size(dat.binned,2));
standard_same_context_err = NaN(size(dat.timebins,2), 2, size(dat.binned,2));

oddball_same_context = NaN(size(dat.timebins,2),2, size(dat.binned,2));
oddball_diff_context = NaN(size(dat.timebins,2), 2, size(dat.binned,2));
oddball_same_context_err = NaN(size(dat.timebins,2), 2, size(dat.binned,2));
dat.durations_old = dat.durations;
dat.durations = [1; 2];

for i = 1:size(dat.binned,2) % for ech channel
    allDurSameCorr= [];
    fprintf(strcat(['channel #' num2str(dat.channel_numbers(i)) '\n']))


    figMatrix = figure('Position', [1 1 700 500], 'visible', 'off');
    sgtitle(figMatrix, {[ 'sub-' lower(dat.indiv) ' -- sess-' dat.sessN ' -- ' dat.taskType], ...
    ['signal-' signalType ' -- param-' dat.processing_params], ...
    ['binWidth-' num2str(dat.binWidthMs) 'ms -- chann#'  num2str(dat.channel_numbers(i)) '-- testretest-corr' ]},'FontSize',10);

    figTCI = figure('Position', [1 1 800 400]) %, 'visible', 'off');
    sgtitle(figTCI, {[ 'sub-' lower(dat.indiv) ' -- sess-' dat.sessN ' -- ' dat.taskType], ...
    ['signal-' signalType], ...
    ['param-' dat.processing_params], ...
    ['binWidth-' num2str(dat.binWidthMs) 'ms -- chann#'  num2str(dat.channel_numbers(i)) '-- TCI' ]},'FontSize',10);

    for j=1:size(dat.durations,1) % for each duration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %STANDARD
        
    
        % calculate same and cross context correlatio at each time point
        for t=1:size(dat.binned(i).chann.standard.order1,4)
            sameR1 = corr(squeeze(dat.binned(i).chann.standard.order1(j,:,:,t))'); %stim(D,A)x rep x pattern 
            sameR2 = corr(squeeze(dat.binned(i).chann.standard.order2(j,:,:,t))');
            
            sameCorrT(1, t) = mean(sameR1(logical(tril(ones(size(sameR1)),-1))), 'all', 'omitnan');
            sameCorrT(2, t) = mean(sameR2(logical(tril(ones(size(sameR2)),-1))), 'all', 'omitnan');
            
            diffR = corr(squeeze(dat.binned(i).chann.standard.order1(j,:,:,t))', squeeze(dat.binned(i).chann.standard.order2(j,:,:,t))');
            diffCorr(1,t) = mean(diffR, 'all', 'omitnan');
        end
        
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO FINISH   
%         % calculate sam context over time
%         sameCorr1 = NaN(size(dat.binned(i).chann.standard.order1,2), ...
%             size(dat.binned(i).chann.standard.order1,1),...
%             size(dat.binned(i).chann.standard.order1,4));
%         
%         sameCorr2 = NaN(size(dat.binned(i).chann.standard.order2,2),...
%             size(dat.binned(i).chann.standard.order2,1),...
%             size(dat.binned(i).chann.standard.order2,4));
%         
%         for s=1:size(dat.binned(i).chann.standard.order1,2)
%             for t=1:size(dat.binned(i).chann.standard.order1,4)
% 
%                 sameCorr1(s,:,t) = corr(squeeze(dat.binned(i).chann.standard.order1(:,s,:,t))');
%                 sameCorr2(s,:,t) = corr(squeeze(dat.binned(i).chann.standard.order2(:,s,:,t))');
%             
%             end
%         end
%         sameCorr = [sameCorr1; sameCorr2];

        
%         % plot matrix for each duration
%         set(0,'CurrentFigure',figMatrix)
%         subplot(size(dat.durations,1),2,j);
%         imagesc(squeeze(mean(sameCorr,1, 'omitnan'))); 
%         colormap(jet);
%         colorbar;
%         caxis([0 1]);
%         title(strcat(['dur-' num2str(floor(dat.durations(j)*1000))  'ms']));
        
        set(0,'CurrentFigure',figTCI)
        subplot(size(dat.durations,1), 2, j);
        plot(dat.timebins, mean(sameCorrT,1), 'color', [.5 .5 .5]); hold on; 
        plot(dat.timebins, diffCorr);
        xline(0, '-b', 'LineWidth', 1);% xline(dat.durations(j), '-b', 'LineWidth', 1);
%         xlim([0 .2]);
        ylim([-.1 max(mean(sameCorrT,1)) + .1]);
        title(strcat(['standard']));
%         ylim([min(diffCorr)*1.3 max(mean(sameCorr,1))*1.5]);
%         yticks([-0.05 0 0.1 0.2]);

%         allDurSameCorr(:,:,j) = [ squeeze(mean(sameCorr,1, 'omitnan'))];

        standard_same_context(:,j,i)=mean(sameCorrT,1);
        standard_diff_context(:,j,i)=diffCorr;
        standard_same_context_err(:,j,i) = (sameCorrT(1,:)/2 - sameCorrT(2,:)/2).^2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ODDBALL
        
    
        % calculate same and cross context correlatio at each time point
        for t=1:size(dat.binned(i).chann.oddball.order1,4)
            sameR1 = corr(squeeze(dat.binned(i).chann.oddball.order1(j,:,:,t))'); %stim(D,A)x rep x pattern 
            sameR2 = corr(squeeze(dat.binned(i).chann.oddball.order2(j,:,:,t))');
            
            sameCorrT(1, t) = mean(sameR1(logical(tril(ones(size(sameR1)),-1))), 'all', 'omitnan');
            sameCorrT(2, t) = mean(sameR2(logical(tril(ones(size(sameR2)),-1))), 'all', 'omitnan');
            
            diffR = corr(squeeze(dat.binned(i).chann.oddball.order1(j,:,:,t))', squeeze(dat.binned(i).chann.oddball.order2(j,:,:,t))');
            diffCorr(1,t) = mean(diffR, 'all', 'omitnan');
        end
        
        
%         % calculate sam context over time
%         sameCorr1 = NaN(size(dat.binned(i).chann.oddball.order1,2), ...
%             size(dat.binned(i).chann.oddball.order1,3),...
%             size(dat.binned(i).chann.oddball.order1,3));
%         
%         sameCorr2 = NaN(size(dat.binned(i).chann.oddball.order2,2),...
%             size(dat.binned(i).chann.oddball.order2,3),...
%             size(dat.binned(i).chann.oddball.order2,3));
%         
%         for s=1:size(dat.binned(i).chann.oddball.order1,2)
% 
%             sameCorr1(s,:,:) = corr(squeeze(dat.binned(i).chann.oddball.order1(:,s,:,:))');
%             sameCorr2(s,:,:) = corr(squeeze(dat.binned(i).chann.oddball.order2(:,s,:,:))');
%             
% 
%         end
%         sameCorr = [sameCorr1; sameCorr2];

        
%         % plot matrix for each duration
%         set(0,'CurrentFigure',figMatrix)
%         subplot(size(dat.durations,1),2,j+2);
%         imagesc(squeeze(mean(sameCorr,1, 'omitnan'))); 
%         colormap(jet);
%         colorbar;
%         caxis([0 1]);
%         title(strcat(['dur-' num2str(floor(dat.durations(j)*1000))  'ms']));
        
        set(0,'CurrentFigure',figTCI)
        subplot(size(dat.durations,1), 2, j+2);
        plot(dat.timebins, mean(sameCorrT,1), 'color', [.5 .5 .5]); hold on; 
        plot(dat.timebins, diffCorr);
        xline(0, '-b', 'LineWidth', 1); %xline(dat.durations(j), '-b', 'LineWidth', 1);
%         xlim([0 .2]);
        ylim([-.1 max(mean(sameCorrT,1)) + .1]);
        title(strcat(['oddball']));
%         ylim([min(diffCorr)*1.3 max(mean(sameCorr,1))*1.5]);
%         yticks([-0.05 0 0.1 0.2]);

%         allDurSameCorr(:,:,j) = [ squeeze(mean(sameCorr,1, 'omitnan'))];

        oddball_same_context(:,j,i)=mean(sameCorrT,1);
        oddball_diff_context(:,j,i)=diffCorr;
        oddball_same_context_err(:,j,i) = (sameCorrT(1,:)/2 - sameCorrT(2,:)/2).^2;
    end
    set(0,'CurrentFigure',figTCI)
    xlabel('Time from segment onset(s)');
    ylabel('same/cross-context corr');
    
    t = strcat([ figureFolder filesep 'sub-' lower(dat.indiv) filesep 'sess-' dat.sessN filesep dat.taskType filesep dat.trialType filesep ... 
                                'TCI' filesep dat.signalTime filesep  signalType filesep dat.processing_params filesep... 
                                'binWidth-' num2str(dat.binWidthMs) 'ms' filesep 'eachDur'  ]);
                            

    if ~isfolder(t)
        mkdir(t)  
    end
    
    saveas(figTCI, strcat([t filesep 'sub-' lower(dat.indiv) '_ses-' dat.sessN '_task-' lower(dat.taskType) '_trial-' lower(dat.trialType) '_binWidth_' num2str(dat.binWidthMs) 'ms_TCI_chann-' num2str(i) '.png']));
    close(figTCI);
    
    set(0,'CurrentFigure',figMatrix)
    
    t = strcat([ figureFolder filesep 'sub-' lower(dat.indiv) filesep 'sess-' dat.sessN filesep dat.taskType filesep dat.trialType filesep ... 
                                'testRetestCorr' filesep dat.signalTime filesep  signalType filesep dat.processing_params filesep... 
                                'binWidth-' num2str(dat.binWidthMs) 'ms' filesep 'eachDur'  ]);

    if ~isfolder(t)
        mkdir(t)  
    end
    
    saveas(figMatrix, strcat([t filesep 'sub-' lower(dat.indiv) '_ses-' dat.sessN '_task-' lower(dat.taskType) '_trial-' lower(dat.trialType) '_binWidth_' num2str(dat.binWidthMs) 'ms_testRetestCorr_chann-' num2str(i) '.png']));
    close(figMatrix);
    

%     figMatrixAll = figure('Position', [1 1 700 500], 'visible', 'off');
% 
%     imagesc(squeeze(mean(allDurSameCorr,3, 'omitnan'))); 
%     colormap(jet);
%     colorbar;
%     caxis([0 1]);
%     title(strcat( ['Average testretest-corr for all durations' ]));
%     subtitle({[ 'sub-' dat.indiv ' -- sess-' dat.sessN ' -- ' dat.taskType], ...
%     ['signal-' signalType ' -- param-' dat.processing_params], ...
%     ['binWidth-' num2str(dat.binWidthMs) 'ms -- chann#'  num2str(dat.channel_numbers(i))]},'FontSize',10);
% 
%     t = strcat([ figureFolder filesep 'sub-' lower(dat.indiv) filesep 'sess-' dat.sessN filesep dat.taskType filesep dat.trialType filesep ... 
%                                 'testRetestCorr' filesep dat.signalTime filesep  signalType filesep dat.processing_params filesep... 
%                                 'binWidth-' num2str(dat.binWidthMs) 'ms' filesep 'allDurs'  ]);
% 
%     if ~isfolder(t)
%         mkdir(t)  
%     end
%     
%     saveas(figMatrixAll, strcat([t filesep 'sub-' lower(dat.indiv) '_ses-' dat.sessN '_task-' lower(dat.taskType) '_trial-' lower(dat.trialType) '_binWidth_' num2str(dat.binWidthMs) 'ms_testRetestCorr_chann-' num2str(i) '.png']));
%     close(figMatrixAll);


end



L.unique_segs = [dat.binned(1).chann(:).durMs];
L.sr = dat.srEphysDownSr;
L.lag_t = dat.timebins;
L.n_total_segs = [dat.binned(1).chann(:).durN]';
L.channels = dat.channel_numbers;
L.param_string = ['boundary-any_lag_win-0-1'];
L.chnames = cellstr(num2str(L.channels'));
L.boundary = 'any';

figure_directory = strcat([source_folder filesep 'analysis' filesep 'figures' filesep ...
    'sub-' lower(dat.indiv) filesep 'sess-' num2str(dat.sessN) filesep dat.taskType filesep dat.trialType filesep ...
    'TCI' filesep 'stim-aligned' filesep signalType filesep dat.processing_params filesep 'binWidth-' num2str(dat.binWidthMs) 'ms' filesep ...
     L.param_string]);
if ~isfolder(figure_directory)
    mkdir(figure_directory)  
end
output_directory = strcat([source_folder filesep 'analysis' filesep 'data-work' filesep ...
    'sub-' lower(dat.indiv) filesep 'sess-' num2str(dat.sessN) filesep dat.taskType filesep dat.trialType filesep ...
    'TCI' filesep 'stim-aligned' filesep signalType filesep dat.processing_params filesep 'binWidth-' num2str(dat.binWidthMs) 'ms' filesep ...
     L.param_string]);
if ~isfolder(output_directory)
    mkdir(output_directory)  
end 
L.figure_directory = figure_directory;
L.output_directory = output_directory;
L.rampwin = 'hann';
L.rampdur = 0.015625 ;
L.same_context = standard_same_context;
L.diff_context = standard_diff_context;
L.same_context_err = standard_same_context_err;
L.indiv = dat.indiv;
L.sessN = dat.sessN;
L.taskType = dat.taskType;
L.trialType = dat.trialType;

save([output_directory filesep 'sub-' lower(dat.indiv) '_ses-' dat.sessN '_task-' lower(dat.taskType) '_trial-' lower(dat.trialType)  '_Lstruct.mat'], 'L');
end