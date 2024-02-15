function [id, rcc_sq, phi_out] = tw_evalPlaneFits(phi, phiPred, movMeanWinSize)
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




%%
id = nan(nTm,nTr);
rcc_sq = nan(nTm,nTr);
phi_out = nan(size(phi,1), nTm, nTr);

wb = waitbar(0, {sprintf('Wave Fit: Trial %d/%d', 1,nTr), '~~~ >> ??? << ~~~'});

for itrial = 1:nTr
    waitbar(itrial/nTr,wb, {sprintf('Wave Fit: Trial %d/%d', itrial,nTr), '~~~ >> ??? << ~~~'});
    for iT = 1:nTm
        aphi = squeeze(circ_mean(phi(:,movAvgIdx(iT),itrial), [], 2));
        gof = sqrt(mean(cos(phiPred'-aphi)).^2 + mean(sin(phiPred'-aphi)).^2); % mean vector length of residuals (-> ie offset is ignored)
        [~,id(iT,itrial)] = max(gof); % get best fit ID
        bestfit = phiPred(id(iT,itrial),:)'; % get predicted phases for best fit
        rcc_sq(iT,itrial) = sum(sin(aphi-circ_mean(aphi)) .* sin(bestfit-circ_mean(bestfit)))...
            ./sqrt(...
            sum(sin(aphi-circ_mean(aphi)).^2) .* sum(sin(bestfit-circ_mean(bestfit)).^2)); % circular correlation for best fit
        rcc_sq(iT,itrial) = rcc_sq(iT,itrial)^2; 
        phi_out(:,iT,itrial) = aphi;
    end
end
delete(wb);

phi_out = phi_out - circ_mean(phi_out, [], 1); % recenter on zero after moving mean

end