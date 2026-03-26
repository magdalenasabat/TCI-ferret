function triggIdx = myGetTrigIdx(triggerChannel, srEphys)

%threshold triggerChannel (over 600)
th = 2000;
triggerSmp = zeros(1, length(triggerChannel));
triggerSmp(1,triggerChannel>th) = 1;

durTriggerMs = 0.01; %it's 5 ms but in 44,1kHz so round up to 6ms in case of some fukcup
durTriggerSmp = durTriggerMs*srEphys; 
durTriggIdx = zeros(1, length(triggerChannel));
triggIdx = [];
lastIdx = 1; 

% timings = find(gt(2:end)>gt(1:(end-1)))+1; %find samples when n+1 is
% larger than n THIS IS SO MUCH BETTER

for i=1:length(triggerSmp)
    if (triggerSmp(i)==1 && i-lastIdx > durTriggerSmp) || (triggerSmp(i)==1 && lastIdx==1) %% i's a trigger sample (==1) and the distance from last trigger sample > duration of a trigger or it's a trigger sample and it's a first trigger sample%         
        durTriggIdx(i) = 1;
        triggIdx = [triggIdx, i];
        lastIdx = i+1;
    end
end
end