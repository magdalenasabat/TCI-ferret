function noise=myCorrReps(binned)

noise = [];
nrep =  size(binned(1).chann(1).order1,2);
for chan=1:length(binned)
    c = [];
    for dur=1:length(binned(chan).chann)
        r1 = corrcoef(binned(chan).chann(dur).order1(:,1:2:nrep,:), ...
                          binned(chan).chann(dur).order1(:,2:2:nrep,:));
        r2 = corrcoef(binned(chan).chann(dur).order2(:,1:2:nrep,:), ...
                          binned(chan).chann(dur).order2(:,2:2:nrep,:));
        c(dur) = sqrt(mean([r1(1,2),r2(1,2)]));
    end
    noise(chan) = mean(c);

end