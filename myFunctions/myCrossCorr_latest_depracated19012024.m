function L = myCrossCorr_latest(dat, source_folder, plot_figure, smoothing, saveDat )
signalType=dat.signalType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% intitialise

if ~exist('plot_figure','var')
  plot_figure = 'False';
end
if ~exist('saveDat','var')
  saveDat = 1;
end

figureFolder = strcat([ source_folder filesep 'analysis' filesep 'figures']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cleanup


fprintf('calculating correlations  \n')

same_context = NaN(size(dat.timebins,2), size(dat.durations,1), size(dat.binned,2));
diff_context = NaN(size(dat.timebins,2), size(dat.durations,1), size(dat.binned,2));
same_context_err = NaN(size(dat.timebins,2), size(dat.durations,1), size(dat.binned,2));


for i = 1:size(dat.binned,2) % for ech channel
    allDurSameCorr= [];
    fprintf(strcat(['channel #' num2str(dat.channel_numbers(i)) '\n']))

% 
%     figMatrix = figure('Position', [1 1 700 500], 'visible', 'off');
%     sgtitle(figMatrix, {[ 'sub-' lower(dat.indiv) ' -- sess-' dat.sessN ' -- ' dat.taskType], ...
%     ['signal-' signalType ], ...
%     ['binWidth-' num2str(dat.binWidthMs) 'ms -- chann#'  num2str(dat.channel_numbers(i)) '-- testretest-corr' ]},'FontSize',10);

    figTCI = figure('Position', [1 1 400 800], 'visible', 'off');
    sgtitle(figTCI, {[ 'sub-' lower(dat.indiv) ' -- sess-' dat.sessN ' -- ' dat.taskType], ...
    ['signal-' signalType], ...
    ['binWidth-' num2str(dat.binWidthMs) 'ms -- chann#'  num2str(dat.channel_numbers(i)) '-- TCI' ]},'FontSize',10);

    for j=1:size(dat.durations,1) % for each duration
        
        sameCorrT = NaN(2, size(dat.binned(i).chann(j).order1,3));
        diffCorr = NaN(1, size(dat.binned(i).chann(j).order1,3));

        % calculate same and cross context correlatio at each time point
        for t=1:size(dat.binned(i).chann(j).order1,3)
            switch signalType
                case 'stimuli'
                    sameCorr1 =  mean(corr(squeeze(dat.binned(i).chann(j).order1)),1);
                    sameCorr2 =  mean(corr(squeeze(dat.binned(i).chann(j).order2)),1);
                otherwise
                
                    sameR1 = corr(dat.binned(i).chann(j).order1(:,:,t));
                    sameR2 = corr(dat.binned(i).chann(j).order2(:,:,t));
            end
            
            sameCorrT(1, t) = mean(sameR1(logical(tril(ones(size(sameR1)),-1))), 'all', 'omitnan');
            sameCorrT(2, t) = mean(sameR2(logical(tril(ones(size(sameR2)),-1))), 'all', 'omitnan');
            
            diffR = corr(dat.binned(i).chann(j).order1(:,:,t), dat.binned(i).chann(j).order2(:,:,t));
            diffCorr(1,t) = mean(diffR, 'all', 'omitnan');
        end
        
        
        % calculate sam context over time
        sameCorr1 = NaN(size(dat.binned(i).chann(j).order1,1),size(dat.binned(i).chann(j).order1,2),size(dat.binned(i).chann(j).order1,2));
        sameCorr2 = NaN(size(dat.binned(i).chann(j).order2,1),size(dat.binned(i).chann(j).order2,2),size(dat.binned(i).chann(j).order2,2));
        
        switch signalType
            case 'stimuli'
                sameCorr1 =  mean(corr(squeeze(dat.binned(i).chann(j).order1)),1);
                sameCorr2 =  mean(corr(squeeze(dat.binned(i).chann(j).order2)),1);
                
            otherwise
                for s=1:size(dat.binned(i).chann(j).order1,1)

                    sameCorr1(s,:,:) = corr(squeeze(dat.binned(i).chann(j).order1(s,:,:))');
                    sameCorr2(s,:,:) = corr(squeeze(dat.binned(i).chann(j).order2(s,:,:))');
            

                end
        end

        sameCorr = [sameCorr1; sameCorr2];

        
%         % plot matrix for each duration
%         set(0,'CurrentFigure',figMatrix)
%         subplot(ceil(size(dat.durations,1)/3),ceil(size(dat.durations,1)/2),j);
%         imagesc(squeeze(mean(sameCorr,1, 'omitnan'))); 
%         colormap(jet);
%         colorbar;
%         caxis([0 1]);
%         title(strcat(['dur-' num2str(floor(dat.durations(j)*1000))  'ms']));
        
        set(0,'CurrentFigure',figTCI)
        subplot(size(dat.durations,1), 1, j);
        plot(dat.timebins, mean(sameCorrT,1), 'color', [.5 .5 .5]); hold on; 
        plot(dat.timebins, diffCorr);
        xline(0, '-b', 'LineWidth', 1); xline(dat.durations(j), '-b', 'LineWidth', 1);
        title(strcat(['dur-' num2str(floor(dat.durations(j)*1000))  'ms']));
%         ylim([min(diffCorr)*1.3 max(mean(sameCorr,1))*1.5]);
%         yticks([-0.05 0 0.1 0.2]);

%         allDurSameCorr(:,:,j) = [ squeeze(mean(sameCorr,1, 'omitnan'))];

        same_context(:,j,i)=mean(sameCorrT,1);
        diff_context(:,j,i)=diffCorr;
        same_context_err(:,j,i) = (sameCorrT(1,:)/2 - sameCorrT(2,:)/2).^2;

    end
    set(0,'CurrentFigure',figTCI)
    xlabel('Time from segment onset(s)');
    ylabel('same/cross-context corr');
    
    t = strcat([ dat.figure_directory filesep 'eachDur'  ]);
                            

    if ~isfolder(t)
        mkdir(t)  
    end
    
    saveas(figTCI, strcat([t filesep 'sub-' lower(dat.indiv) '_ses-' dat.sessN '_task-' lower(dat.taskType) '_trial-' lower(dat.trialType) '_binWidth_' num2str(dat.binWidthMs) 'ms_TCI_chann-' num2str(i) '.png']));
    close(figTCI);
    
% %     set(0,'CurrentFigure',figMatrix)
%     
%     t = strcat([ figureFolder filesep 'sub-' lower(dat.indiv) filesep 'sess-' dat.sessN filesep dat.taskType filesep dat.trialType filesep ... 
%                                 'testRetestCorr' filesep dat.signalTime filesep  signalType filesep... 
%                                 'binWidth-' num2str(dat.binWidthMs) 'ms' filesep 'eachDur'  ]);
% 
%     if ~isfolder(t)
%         mkdir(t)  
%     end
%     
%     saveas(figMatrix, strcat([t filesep 'sub-' lower(dat.indiv) '_ses-' dat.sessN '_task-' lower(dat.taskType) '_trial-' lower(dat.trialType) '_binWidth_' num2str(dat.binWidthMs) 'ms_testRetestCorr_chann-' num2str(i) '.png']));
%     close(figMatrix);
    

%     figMatrixAll = figure('Position', [1 1 700 500], 'visible', 'off');
% 
%     imagesc(squeeze(mean(allDurSameCorr,3, 'omitnan'))); 
%     colormap(jet);
%     colorbar;
%     caxis([0 1]);
%     title(strcat( ['Average testretest-corr for all durations' ]));
%     subtitle({[ 'sub-' dat.indiv ' -- sess-' dat.sessN ' -- ' dat.taskType], ...
%     ['signal-' signalType], ...
%     ['binWidth-' num2str(dat.binWidthMs) 'ms -- chann#'  num2str(dat.channel_numbers(i))]},'FontSize',10);
% 
%     t = strcat([ figureFolder filesep 'sub-' lower(dat.indiv) filesep 'sess-' dat.sessN filesep dat.taskType filesep dat.trialType filesep ... 
%                                 'testRetestCorr' filesep dat.signalTime filesep  signalType filesep... 
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
L.sr = dat.sr;
L.lag_t = dat.timebins;
L.n_total_segs = [dat.binned(1).chann(:).durN]';
L.channels = dat.channel_numbers;
L.param_string = ['boundary-any_lag_win-0-1'];
L.chnames = cellstr(num2str(L.channels'));
L.boundary = 'any';



L.rampwin = 'hann';
L.rampdur = 0.015625 ;
L.same_context = same_context;
L.diff_context = diff_context;
L.same_context_err = same_context_err;
L.indiv = dat.indiv;
L.sessN = dat.sessN;
L.taskType = dat.taskType;
L.trialType = dat.trialType;

if saveDat

    L.figure_directory =      strcat([  dat.figure_directory filesep L.param_string]);
    L.output_directory =      strcat([  dat.output_directory filesep L.param_string]); 
    
    if ~isfolder(L.figure_directory)
        mkdir(L.figure_directory)  
    end
    if ~isfolder(L.output_directory)
        mkdir(L.output_directory)  
    end        

    save([L.output_directory filesep 'sub-' lower(dat.indiv) '_ses-' dat.sessN '_task-' lower(dat.taskType) '_trial-' lower(dat.trialType)  '_Lstruct.mat'], 'L', '-v7.3');
    
end
end