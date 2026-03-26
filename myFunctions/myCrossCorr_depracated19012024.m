function L = myCrossCorr(d, save_folder, plot_figure, smoothing)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% intitialise

if ~exist('plot_figure','var')
  plot_figure = 'False';
end

stimuli_names = d.stimuli_names;
indiv = d.indiv;
channel_numbers = d.channel_numbers;
srEphys = d.srEphys;
spikes = d.spikes;
triggers= d.triggers;
% lfp = d.lfp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cleanup

if strcmp(stimuli_names{1}, 'silence - buffer init')
    stimuli_names(1) = [];
end
% if strcmpi(dataMat.trialType, 'Active')
%     stimuli_names(find(x==find(strcmpi(y,dataMat.targetName)))) = [];
% %     stimuliFileTCI(find(x==find(strcmpi(y,dataMat.targetName)))) = [];
%     triggers(find(x==find(strcmpi(y,dataMat.targetName)))) = [];
% end
% in = find(strcmpi(stimuli_names,'speech-french_500ms_part14_31.25ms_part5.wav') | strcmpi(stimuli_names,dataMat.targetName));
% if ~isempty(in)
%     stimuli_names(in) = [];
% %     stimuliFileTCI(in) = [];
%     triggers(in) = [];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% create strict with triggers for each stimuli
fprintf('getting triggers  \n')
%remove target presentations
% [x,y]  = grp2idx(stimuli_names);


durStimuli = myGetDuration(stimuli_names);
durStimuli = cellfun(@(x) str2num(x)/1000, durStimuli);
[durIdx,~, dur]  = grp2idx(durStimuli);
[stimIdx,namesTCI] = grp2idx(stimuli_names);

pattern1 = stimIdx(1:length(namesTCI));
pattern2 = stimIdx(length(namesTCI)+1:(length(namesTCI)+length(namesTCI)));

pattern1onsets = strfind(stimIdx',pattern1');
pattern2onsets = strfind(stimIdx',pattern2');

sz = min(length(pattern1onsets),length(pattern2onsets));

pattern1onsets=pattern1onsets(1:sz);
pattern2onsets=pattern2onsets(1:sz);

A = pattern1onsets+pattern1-1;
A = reshape(A,[],1) ;

B = pattern2onsets+pattern1-1;
B = reshape(B,[],1) ;

%%
fprintf('calculating correlations  \n')
beforeT= 0.1*srEphys;
afterT= 0.5*srEphys; 
dt = .001*srEphys;
timebins = -beforeT:dt:afterT;
% time = (-beforeT:dt:afterT)/srEphys;
t1 = tic;
same_context = NaN(size(timebins,2), size(dur,1), size(unique(d.spikes(:,2)),1));
diff_context = NaN(size(timebins,2), size(dur,1), size(unique(d.spikes(:,2)),1));
same_context_err = NaN(size(timebins,2), size(dur,1),size(unique(d.spikes(:,2)),1));

%for tests
% same_context = NaN(afterT/10, size(dur,1), size(unique(d.spikes(:,2)),1));
% diff_context = NaN(afterT/10, size(dur,1), size(unique(d.spikes(:,2)),1));
% same_context_err = NaN(afterT/10, size(dur,1),size(unique(d.spikes(:,2)),1));
time = (1:1:afterT/10)/(srEphys/10);

%for each channel
for i = 1:size(unique(d.spikes(:,2)))
    fprintf(strcat(['channel #' num2str(channel_numbers(i)) '\n']))
    XX = spikes(spikes(:,2)==i,1); 
%     
%     spikesIdx = find(data> (mean(data)+ std(data)*3) ...
%         | (data < (mean(data)- std(data)*3)));
%     
%     XX = spikesIdx; %for each duration
%     XX = lfp(i,:);
    if plot_figure == true; figure('Position', [1 1 400 800], 'visible', 'off');sgtitle(['chann #' num2str(channel_numbers(i))]); end
    
    for j=1:size(dur,1)        
        order1 = [];
        order2 = [];
        order1Idx = [];
        order2Idx = [];
        
        indexesDur = find(durIdx == j); %get indexes of all stim of a given dur
        
        x1 = length(intersect(pattern1,indexesDur));
        x2 = length(intersect(pattern2,indexesDur));
        
        id1=intersect(indexesDur, A);
        id2=intersect(indexesDur, B);
        
        order1Idx = triggers(id1);
        order1Idx = reshape(order1Idx, [x1, length(order1Idx)/x1]);
        
        order2Idx = triggers(id2);
        order2Idx = reshape(order2Idx, [x2, length(order2Idx)/x2]);
      
        order1 = NaN(size(order1Idx,1), size(order1Idx,2) , length(timebins));
        order2 = NaN(size(order2Idx,1), size(order2Idx,2) , length(timebins));

%for tests
% %         afterT = floor(dur(j)*44100);
%         order1 = NaN(size(order1Idx,1), size(order1Idx,2) , afterT/10);
%         order2 = NaN(size(order2Idx,1), size(order2Idx,2) , afterT/10);
        
        v = stimuli_names(id1);
        v(2,:) = stimuli_names(id2);
        v = v(:,1:x1);
        v=[v(1,:) v(2,:)];
        cool = grp2idx(v);
        cool2 = cool(x1+1:end);
        back=order2Idx;
        
        for k=1:size(order2Idx,2)
            for m=1:size(cool2,1)
                order2Idx(m,k) = back(find(cool2==m),k);
            end         
        end
        
        for k=1:size(order1Idx,1)
            for m=1:size(order1Idx,2)          
%                 order1(k,m,:) = hist(XX(ismember(XX, order1Idx(k,m)-beforeT:1:order1Idx(k,m)+afterT-1))- (order1Idx(k,m)-beforeT) , timebins);
%                 order2(k,m,:) = hist(XX(ismember(XX, order2Idx(k,m)-beforeT:1:order2Idx(k,m)+afterT-1))- (order2Idx(k,m)-beforeT) , timebins);
                if smoothing
                    order1(k,m,:) = smoothdata( ...
                        hist(XX(XX > order1Idx(k,m) - beforeT & XX < order1Idx(k,m) + afterT) - order1Idx(k,m) , timebins), ...
                        'gaussian',5);
                    order2(k,m,:) = smoothdata( ...
                        hist(XX(XX > order2Idx(k,m) - beforeT & XX < order2Idx(k,m) + afterT) - order2Idx(k,m) , timebins), ...
                        'gaussian',5);
                else
                    order1(k,m,:) = hist(XX(XX > order1Idx(k,m) - beforeT & XX < order1Idx(k,m) + afterT) - order1Idx(k,m) , timebins);
                    order2(k,m,:) = hist(XX(XX > order2Idx(k,m) - beforeT & XX < order2Idx(k,m) + afterT) - order2Idx(k,m) , timebins);
%for tests
%     
%                     order1(k,m,:) = downsample(XX(order1Idx(k,m): order1Idx(k,m) + afterT-1),10);
%                     order2(k,m,:) = downsample(XX(order2Idx(k,m): order2Idx(k,m) + afterT-1),10);
                end
            end
        end
%         figure();plot(mean(squeeze(order1(:,1,:)),1))        
        % corr same context
        sameCorr = NaN(2,size(order1,3));
        diffCorr = NaN(1,size(order1,3));
        
        for t=1:size(order1,3)
            sameR1 = corr(order1(:,:,t));
            sameR2 = corr(order2(:,:,t));
%             if dur(j) > .45
%                 figure();
% %                 title();
%                 subplot(1,2,1)
%                 imagesc(sameR1);
%                 subplot(1,2,2)
%                 imagesc(sameR2);
%                 colorbar;
% %                 title(num2str(dur(j)));
%             end

            sameCorr(1, t) = mean(nonzeros(tril(sameR1,-1)), 'all');
            sameCorr(2, t) = mean(nonzeros(tril(sameR2,-1)), 'all');
            
            
            diffR = corr(order1(:,:,t), order2(:,:,t));
            diffCorr(1,t) = mean(diffR, 'all');

        end 

            
        same_context(:,j,i)=mean(sameCorr,1);
        diff_context(:,j,i)=diffCorr;
        same_context_err(:,j,i) = (sameCorr(1,:)/2 - sameCorr(2,:)/2).^2;;
        time=timebins/srEphys;
        if plot_figure == true
            subplot(size(dur,1), 1, j);
            plot(time, mean(sameCorr,1)); hold on; plot(time, diffCorr);
            xline(0, '-r', 'LineWidth', 1); xline(dur(j), '-r', 'LineWidth', 1);
            title([num2str(dur(j))]); 
            xlim([-beforeT/srEphys afterT/srEphys]);
            ylim([-0.05 max(max(mean(sameCorr,1))*1.3, 0.05)]);
            yticks([-0.05 0 0.1 0.2]);
        end
    end
    if plot_figure == true
        xlabel('Time from segment onset');
        ylabel('noise ceil / cross-context correlation ');
        if smoothing
            saveas(gcf,[save_folder '_TCI-'  num2str(channel_numbers(i)) '_smoothed.png']);
        else
            saveas(gcf,[save_folder '_TCI-'  num2str(channel_numbers(i)) '.png']);
            
%         close;
    end

end

% time = (-beforeT:dt:afterT)/srEphys;
loc = save_folder ;
tbl = tabulate(grp2idx(durStimuli(1:length(pattern1))));

L.unique_segs = dur;
L.sr = 1000.;
L.lag_t = time'; %????????
L.n_total_segs = [tbl(:,2)];
L.channels = 1:size(same_context, 3);
L.param_string = ['boundary-any_lag_win-0-1'];
L.chnames = cellstr(num2str(L.channels'));
L.boundary = 'any';
L.figure_directory = [loc 'figures/v1/boundary-any_lag_win-0-1' ];
L.output_directory = [loc 'results/v1/boundary-any_lag_win-0-1'];
L.rampwin = 'hann';
L.rampdur = 0.015625 ;
L.same_context = same_context;
L.diff_context = diff_context;
L.same_context_err = same_context_err;
L.indiv = indiv;

% save([loc 'L.mat'], 'L');
end