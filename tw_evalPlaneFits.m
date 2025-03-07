function [id, rcc_sq, rcc_sq_rand, phi_out] = tw_evalPlaneFits(phi, phiPred, movMeanWinSize, nIter, shuffleTrials, shTm)
% compares predicted phases in phiPred to actual phases in phi, after
% applying moving-mean of window size movMeanWinSize (in samples, not ms);
% phases in phi must be relative (e.g. to mean phase across electrode) if
% moving mean is requested
% Returned id is the index of the best fit for each time-point, with the same
% [nTime x nTrial] dimensions as input

%%
nTr = size(phi,3);
nTm = size(phi,2);

%%
if nargin < 3 || isempty(movMeanWinSize)
    movAvgIdx = @(t) t;
else
    movAvgIdx = @(t) (1:nTm) > t-movMeanWinSize/2 & (1:nTm) < t+movMeanWinSize/2;
end
if nargin < 4 || isempty(nIter)
    nIter = 0;
end

%% Random shuffling:
if size(phi,1) >= 10
    for iter = 1:nIter
        phi_rand(:,:,:,iter) = phi(randperm(size(phi,1)),:,:);
    end
else
    % for small channel counts (mostly modeling) we make sure that no two shufflings are identical
    assert(nIter < factorial(size(phi,1)),...
        'Number of iterations is greater than number of possible permutations')
    rperm = perms(1:size(phi,1));
    rperm = rperm(randsample(1:size(rperm,1), nIter),:);
    
    for iter = 1:nIter
        phi_rand(:,:,:,iter) = phi(rperm(iter,:),:,:);
    end
end

phi_rand = phi_rand(:,:,shuffleTrials,:); % limit to the trialset we want

%%
id = nan(nTm,nTr);
rcc_sq = nan(nTm,nTr);
rcc_sq_rand = nan(nTm,numel(shuffleTrials),nIter);
phi_out = nan(size(phi,1), nTm, nTr);

wb = waitbar(0, {sprintf('Wave Fit: Timepoint %d/%d', 0,nTm), '~~~ >> ??? << ~~~'});

phiPred = repmat(phiPred, 1,1,nTr); % expand predicted phases to trial dimension
phiPred = permute(phiPred, [2 3 1]); % move fit dimension to the end

for iT = 1:nTm
    waitbar(iT/nTm,wb, {sprintf('Wave Fit: Timepoint %d/%d', iT, nTm), '~~~ >> ??? << ~~~'});
    
    aphi = squeeze(mean(phi(:,movAvgIdx(iT),:), 2));
    aphi = aphi./abs(aphi); % discard magnitude of moving mean angle
    
    [id(iT,:), rcc_sq(iT,:)] = singleFitEval(aphi, phiPred);
    phi_out(:,iT,:) = aphi;
end

for iT = shTm
    waitbar(iT/nTm, wb, {sprintf('Estimating random distribution: Timepoint %d/%d', iT, nTm), '~~~ >> ??? << ~~~'});

    for iter = 1:nIter
                
        rphi = squeeze(mean(phi_rand(:,movAvgIdx(iT),:,iter), 2));
        rphi = rphi./abs(rphi); % discard magnitude of moving mean angle
        
        [~, rcc_sq_rand(iT,:,iter)] = singleFitEval(rphi, phiPred);
    end
end

    

delete(wb);

phi_out = angle(phi_out./mean(phi_out,1)); % recenter on zero after moving mean

end
%%
function [id, rcc_sq] = singleFitEval(aphi, phiPred)
gof = squeeze(abs(mean(phiPred./aphi)));
if iscolumn(gof) % happens if nTr = 1
    gof = gof';
end
[~,id] = max(gof,[],2); % get best fit ID

bestfit = nan(size(aphi));
for iTr = 1:numel(id)
    bestfit(:,iTr) = squeeze(phiPred(:,iTr,id(iTr))); % get predicted phases for best fit 
end


rcc_sq = sum(cSinVar(aphi) .* cSinVar(bestfit))...
    ./sqrt(...
    sum(cSinVar(aphi).^2) .* sum(cSinVar(bestfit).^2)); % circular correlation for best fit

rcc_sq = rcc_sq.^2;
end
%%
function out = cSinVar(in)
% computes the expression sin(angle(in) - circ_mean(angle(in))) for complex
% input (must all have unit magnitude)
out = in./mean(in);
out = imag(out./abs(out));

end