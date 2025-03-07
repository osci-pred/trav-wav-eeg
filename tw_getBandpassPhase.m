function [phi, pow, eegFilt] = tw_getBandpassPhase(raw, sr, freq)
% bandpass filter raw data and return phases
% raw is raw eeg signal as [nChan x nTime x nTrials]
% freq is band as [min max] in Hz

%%
nTr = size(raw,3);
nTm = size(raw,2);
nElec = size(raw,1);

%% Filter the data
bpFilt = designfilt('bandpassfir', 'StopbandFrequency1', 0.85*freq(1), 'PassbandFrequency1', freq(1),...
    'PassbandFrequency2', freq(2), 'StopbandFrequency2', 1.25*freq(2),...
    'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', sr);


eegFilt = nan(nElec, nTm, nTr);

try
    for itrial = 1:nTr
        for iElec = 1:nElec
            eegFilt(iElec,:,itrial) = filtfilt(bpFilt, raw(iElec,:,itrial));
        end
    end
catch errmsg
    if strcmpi(errmsg.identifier, 'signal:filtfilt:InvalidDimensionsDataShortForFiltOrder')
        rawPad = cat(2,raw,raw,raw);
        eegFilt = nan(nElec, 3*nTm, nTr);
        for itrial = 1:nTr
            for iElec = 1:nElec
                eegFilt(iElec,:,itrial) = filtfilt(bpFilt, rawPad(iElec,:,itrial));
            end
        end
        eegFilt = eegFilt(:,nTm+1:2*nTm,:);
    else
        throw(errmsg)
    end
end


%% extract phase:
for itrial = 1:nTr
    h = hilbert(eegFilt(:,:,itrial)');
    phi(:,:,itrial) = angle(h)';
    pow(:,itrial) = mean(abs(h),2);
end

phi = phi - circ_mean(phi, [], 1); % subtract mean across elecs per timepoint
% -> now we can average over time (within short windows)


phi = exp(1i.*phi);

end