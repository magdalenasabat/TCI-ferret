function L = myCrossCorrOLD(dataMat, dataEphys,srEphys, triggIdx, plot_figure )

if ~exist('plot_figure','var')
  plot_figure = 'False';
end

fprintf('getting mean \n')
globalMean = mean(dataEphys,1);

% pring filtering
%% filter 
fprintf('filtering ephys (that might take a second) \n')
% filter 300-8000
f1 = 300; f2 = 6000;
f1 = f1/srEphys*2;
f2 = f2/srEphys*2;
[b,a] = ellip(4,.5,20,[f1 f2]);
filterParams = [b;a];

% remove the mean
x = double(dataEphys);
x = x - globalMean;
x=x';

dataEphysFiltered =filtfilt(filterParams(1,:),filterParams(2,:),x);
dataEphysFiltered = dataEphysFiltered';

%% create strict with triggers for each stimuli
fprintf('getting triggers  \n')
%remove target presentations
[x,y]  = grp2idx(dataMat.stimuliName);

stimuliNameTCI = dataMat.stimuliName;
stimuliFileTCI = dataMat.stimuliFile;
triggIdxTCI = triggIdx;

if strcmpi(dataMat.trialType, 'Active')
    stimuliNameTCI(find(x==find(strcmpi(y,dataMat.targetName)))) = [];
    stimuliFileTCI(find(x==find(strcmpi(y,dataMat.targetName)))) = [];
    triggIdxTCI(find(x==find(strcmpi(y,dataMat.targetName)))) = [];
end
in = find(strcmpi(stimuliNameTCI,'speech-french_500ms_part14_31.25ms_part5.wav') | strcmpi(stimuliNameTCI,dataMat.targetName));
if ~isempty(in)
    stimuliNameTCI(in) = [];
    stimuliFileTCI(in) = [];
    triggIdxTCI(in) = [];
end


durStimuli = cellfun(@length, stimuliFileTCI);
[durIdx,~, durTCI]  = grp2idx(durStimuli);
[stimIdx,namesTCI] = grp2idx(stimuliNameTCI);

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
fprintf('getting spikes  \n')
beforeT= 0.0*srEphys;
afterT= 0.4*srEphys; 
dt = .001*srEphys;
timebins = 0:dt:(beforeT + afterT);
time = (-beforeT:dt:afterT)/srEphys;
t1 = tic;
same_context = NaN(size(timebins,2), size(durTCI,1), size(dataEphysFiltered,1));
diff_context = NaN(size(timebins,2), size(durTCI,1), size(dataEphysFiltered,1));
same_context_err = NaN(size(timebins,2), size(durTCI,1), size(dataEphysFiltered,1));
%for each channel
for i=1:size(dataEphysFiltered,1)
    tic
    fprintf(['channel #' num2str(i) '\n'])
    data = dataEphysFiltered(i,:);
%     
    spikesIdx = find(data> (mean(data)+ std(data)*3) ...
        | (data < (mean(data)- std(data)*3)));
    
    XX = spikesIdx; %for each duration
    
    if strcmp(plot_figure,'True'); figure('Position', [1 1 400 800]);sgtitle(['chann #' num2str(nThisChannels(i))]); end
    
    for j=1:size(durTCI,1)        
        order1 = [];
        order2 = [];
        order1Idx = [];
        order2Idx = [];
        
        indexesDur = find(durIdx ==j); %get indexes of all stim of a given dur
        
        x1 = length(intersect(pattern1,indexesDur));
        x2 = length(intersect(pattern2,indexesDur));
        
        id1=intersect(indexesDur, A);
        id2=intersect(indexesDur, B);
        
        order1Idx = triggIdxTCI(id1);
        order1Idx = reshape(order1Idx, [x1, length(order1Idx)/x1]);
        
        order2Idx = triggIdxTCI(id2);
        order2Idx = reshape(order2Idx, [x2, length(order2Idx)/x2]);
      
        order1 = NaN(size(order1Idx,1), size(order1Idx,2) , length(timebins));
        order2 = NaN(size(order2Idx,1), size(order2Idx,2) , length(timebins));
        
        v = stimuliNameTCI(id1);
        v(2,:) = stimuliNameTCI(id2);
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
                order1(k,m,:) = hist(XX(ismember(XX, order1Idx(k,m)-beforeT:1:order1Idx(k,m)+afterT-1))- (order1Idx(k,m)-beforeT) , timebins);
                order2(k,m,:) = hist(XX(ismember(XX, order2Idx(k,m)-beforeT:1:order2Idx(k,m)+afterT-1))- (order2Idx(k,m)-beforeT) , timebins);
            end
        end
                
        % corr same context
        sameCorr = NaN(2,size(order1,3));
        diffCorr = NaN(1,size(order1,3));
        
        for t=1:size(order1,3)
            sameR1 = corr(order1(:,:,t));
            sameR2 = corr(order2(:,:,t));
            
            sameCorr(1, t) = mean(nonzeros(tril(sameR1,-1)), 'all');
            sameCorr(2, t) = mean(nonzeros(tril(sameR2,-1)), 'all');
            
            diffR = corr(order1(:,:,t), order2(:,:,t));
            diffCorr(1,t) = mean(diffR, 'all');
            
        end 
        same_context(:,j,i)=mean(sameCorr,1);
        diff_context(:,j,i)=diffCorr;
        same_context_err(:,j,i) = (sameCorr(1,:)/2 - sameCorr(2,:)/2).^2;;
        
        if strcmp(plot_figure,'True')
            subplot(size(durTCI,1), 1, j);
            plot(time, mean(sameCorr,1)); hold on; plot(time, diffCorr);
            xline(0, '-r', 'LineWidth', 1); xline(durTCI(j)/44100, '-r', 'LineWidth', 1);
            title([num2str(durTCI(j)/44100)]); 
            xlim([-beforeT/srEphys afterT/srEphys]);
            ylim([-0.05 max(max(mean(sameCorr,1))*1.3, 0.05)]);
            yticks([-0.05 0 0.1 0.2]);
        end
    end
    if strcmp(plot_figure,'True')
        xlabel('Time from segment onset');
        ylabel('noise ceil / cross-context correlation ');
        saveas(gcf,[mainLocation 'ANALYSIS/myTCI/figures/corr/chann' num2str(nThisChannels(i)) '_TCI.png']);
    end
    close;
    toc
end

time = (-beforeT:dt:afterT)/srEphys;
loc = ['/home' filesep 'magdalena' filesep 'Documents' filesep 'CST' filesep 'CST_PC' filesep 'ANALYSIS' filesep 'myTCI' filesep] ;
tbl = tabulate(grp2idx(durStimuli(1:length(pattern1))));

L.unique_segs = [15.625 31.25 62.5 125 250 500];
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

% save([loc 'L.mat'], 'L');
end