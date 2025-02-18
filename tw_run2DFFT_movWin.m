function out = tw_run2DFFT_movWin(raw, time, varargin)
% run the 2DFFT analysis using a moving window.
% input 'raw' should be EEG as [nChan x nTime x nTrials];
% time should be time vector in seconds

p = inputParser();
p.addParameter('WindowSize', 500); % in ms
p.addParameter('TimeStep', 100); % in ms

p.addParameter('NumIterShuffle', 10);

p.parse(varargin{:});

nChan = size(raw,1);

raw = double(raw);
sr = round(1/diff(time(1:2)));

wsz = round(sr.*(p.Results.WindowSize)./1e3);
tmStep = round(sr.*(p.Results.TimeStep)./1e3);

[fwMax, bwMax, freq, stepTime] = run_main(raw, sr, wsz, tmStep);

[fwMaxBoot, bwMaxBoot] = deal(nan([size(fwMax) p.Results.NumIterShuffle]));

wb = waitbar(0, 'Running random shuffling iterations...');
for iBoot = 1:p.Results.NumIterShuffle
    waitbar(iBoot/p.Results.NumIterShuffle, wb, 'Running random shuffling iterations...');
    rawShuffled = raw(randperm(nChan),:,:);
    [fwMaxBoot(:,:,:,iBoot), bwMaxBoot(:,:,:,iBoot)] = run_main(rawShuffled, sr, wsz, tmStep);
end
delete(wb);

fwMaxBoot = mean(fwMaxBoot,4);
bwMaxBoot = mean(bwMaxBoot,4);
% midlineBoot = mean(midlineBoot,4);

out.fwRaw = fwMax;
out.bwRaw = bwMax;
out.fw = 10.*log10(fwMax./fwMaxBoot);
out.bw = 10.*log10(bwMax./bwMaxBoot);
out.freq = freq;
out.t = time(stepTime);

end
%%
function [fwMax, bwMax, freq, tmVec] = run_main(in, sr, wsz, tmStep)

[~, nTm, nTr] = size(in);

[tmp,~,~,freq] = tw_wavesHunterAllFreqs(in(:,1:wsz,1),sr);

tmVec = wsz/2+1:tmStep:nTm-wsz/2;

[fwMax, bwMax, midline] = deal(nan(numel(tmp), numel(tmVec), nTr));

for iTr = 1:nTr
    
    for iTm = 1:numel(tmVec)
        
        rawTmIdx = tmVec(iTm);
        
        tWin = rawTmIdx + [-0.5 0.5].*wsz;
        
        [fwMax(:,iTm,iTr),bwMax(:,iTm,iTr),midline(:,iTm,iTr)] = tw_wavesHunterAllFreqs(in(:,tWin(1):tWin(2),iTr),sr);
        
    end
    
end
end