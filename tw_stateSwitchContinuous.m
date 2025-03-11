function [state, rTrans, param] = tw_stateSwitchContinuous(wav, varargin)
% Continuous analysis of state switches (not pre-defining FW and BW states)
% 
% input:
%   wav - output of tw_fitPlaneEEG; contains wave fits
%
% optional input (as named arg list):
%
% 'TransThresh'             threshold (in rad phase-diff per ms) to detect
%                           state transitions.
%                           !! This will depend on the parameters used in
%                           the wave fitting before (in particular window
%                           length & frequency). Adjust through visual
%                           check of the data !!
% 'TransMinDist'            min distance (in rad) to be covered by a
%                           state-transition.
% 'TransMaxDur'             max duration of the state transition (in ms).
% 'CircCohThresh'           initial threshold for circular coherence in the
%                           classification of states. This is evaluated from 
%                           mean direction of all significant fits inside 
%                           the time-window, so it takes values within [0
%                           1].
% 'StateMinDur'             min duration of state (in ms).
% 'StateInterpMaxDur'       max duration of gap (non-significant fit) to
%                           interpolate over during stable state.
% 'ShiftToTransMaxTime'     max distance (in ms) to shift state onset/offset to
%                           closest rapid transition (after classification)
% 'DownSampleRate'          rate in Hz to down-sample to in order to speed
%                           up classification process.


p = inputParser();
p.addParameter('TransThresh', 0.2);
p.addParameter('TransMaxDur', 100);
p.addParameter('TransMinDist', pi/4);
p.addParameter('CircCohThresh', 0.85);
p.addParameter('StateMinDur', 100);
p.addParameter('StateInterpMaxDur', 200);
p.addParameter('ShiftToTransMaxTime', 100); 
p.addParameter('DownSampleRate', 50); 

p.parse(varargin{:});


%% Step 1: detect rapid transitions in the wave direction

tStep = diff(wav.t(1:2));
sr = 1/tStep;
nTr = size(wav.wavDir,2);

maxVelEvDur = p.Results.TransMaxDur*1e-3 * sr; % max high-vel event dur in samples

vel = nan(size(wav.wavDir));
rTrans = cell(1,nTr);

for iTr = 1:nTr
    % get velocity as absolute difference between time points, wrapped to
    % [-pi pi]. Values of pi now correspond to jumps of 180 degrees.
    % Scale is radians per ms (we use this only to generalize across
    % sampling rates, the influence of oscillatory frequency and window
    % length is ignored here and needs to be taken into account individually).
    
    vel(:,iTr) = wrapToPi([0; diff(wav.wavDir(:,iTr))])./(tStep*1e3);
    
    ev = bwlabel(abs(vel(:,iTr))' > p.Results.TransThresh); % label high-vel events
    critCheck = @(evID)...
        sum(ev == evID) <= maxVelEvDur... % check for maximum duration of event
        & abs(sum(vel(ev == evID,iTr))*tStep*1e3) >= p.Results.TransMinDist; % check for covered distance
    
    inclEv = arrayfun(critCheck, unique(ev(ev>0)));
    
    % we take the time of the maximum velocity as the event time:
    rTrans{iTr} = arrayfun(@(evId) min(wav.t(abs(vel(:,iTr))'.*(ev == evId) == max(abs(vel(:,iTr))' .* (ev == evId)))), find(inclEv));
    rTrans{iTr} = rTrans{iTr} - 0.5*tStep; % account for left-side diff above
    
%         % plot for sanity check:
%         plot(wav.t, wav.wavDir(:,iTr));
%         hold on
%         for i = find(inclEv)
%             tst = wav.t(find(ev == i, 1, 'first')-1);
%             tend = wav.t(find(ev == i,1, 'last')+1);
%             patch([tst tend tend tst], [-1 -1 1 1].*pi, [1 1 1 1], 'EdgeColor', 'none');
%             alpha 0.3
%         end
%         xline(rTrans{iTr}, 'k:')
end


%% Step 2: detect stable states
% wav.sig is boolean indicating whether the fit was considered significant

srLo = p.Results.DownSampleRate;
[~, tmLo] = resample(wav.wavDir(:,1)', wav.t, srLo); % get time-vector for downsampled signal

maxlen = numel(tmLo);
minlen = round(p.Results.StateMinDur*1e-3 * srLo); 


[xx, yy] = ndgrid(1:numel(tmLo), minlen:maxlen);    


thresh = p.Results.CircCohThresh; % the initial threshold for circular consistency
splitStateThresh = 0.85; % threshold for mean boolean inside the triangular mask 
threshIncr = 0.01; % increment thresh for mask by this amount when splitting states
pad = 0.015; % pad (secs) inward on either side. 
% Pad is only applied to the initial classification, on and offsets
% still adjusted after that based on rapid transitions and state overlaps


interpDur = p.Results.StateInterpMaxDur*1e-3 * srLo; % in samples, how long gap can be to be interpolated

shift2TransMaxTime = p.Results.ShiftToTransMaxTime*1e-3;
state2TransPad = pad; % after moving onset/offset to nearest transition, keep this distance 

wb = waitbar(0, 'State Classification...');

for iTr = 1:nTr
    
    waitbar(iTr/nTr, wb, 'State Classification...');
    
    [wdLo, tmLo] = resample(wav.wavDir(:,iTr)', wav.t, srLo);
    sigLo = interp1(wav.t, double(wav.sig(:,iTr)'), tmLo, 'nearest');
    
    wd = sigLo .* exp(1i.*wdLo);
    % convert to complex, using mag to filter out unsignificant fits
    
    for len = maxlen:-1:minlen
        clust(:,len-minlen+1) = conv(wd, ones(1,len), 'same')./conv(abs(wd), ones(1,len), 'same');
    end
    % square matrix, second dimension is state duration, first is original
    % time
    
    clustMask = abs(clust) > thresh;
    threshMap = thresh + zeros(size(clust)); 
    % threshMap is a map of the thresh to apply for each value of the mask. If
    % we split states somewhere, we differentially adjust values 

    pk = find(any(clustMask, 1), 1, 'last'); % this is the time-scale we start with
    kState = 0;
    
    while pk >= minlen % exit by leaving last step
%         fprintf('\nT = %d', pk)
        clustMask = abs(clust) > threshMap; % update mask with current thresholds
        
        lbl = bwlabel(clustMask(:,pk));
        
        if ~any(lbl)
            pk = pk - 1; % move on directly to the next time-scale
            continue
        end
        
        mag = abs(clust(:, pk));
        ctr = arrayfun(@(x) find(mag.*(lbl==x) == max(mag.*(lbl==x)), 1), unique(lbl(lbl > 0)));
        
        incr = true; % move on after this by default
        for iCtr = 1:numel(ctr)
            
            % create the triangular mask. This represents the state and all
            % its possible subdivisions starting from the timescale T = pk
            dt = ctr(iCtr) + round(pk.*[-0.5 0.5]);
            stateMask1D = (xx >= dt(1) & xx <= dt(2)); % time-only
            stateMask = stateMask1D & (yy <= pk);
            stateMask = stateMask & (abs(yy-pk) >= 2*abs(xx-ctr(iCtr)));
            
            % check that the state doesn't split up into smaller states. This would
            % show as larger gaps in the circular consistency in between this
            % and lower time scales. We apply an arbitrary threshold for the
            % percentage of values above the initial threshold inside the mask.
            splitState = mean(clustMask(stateMask)) < splitStateThresh;
            
            % to-do: add separation criterion by checking outside areas?
            
            if splitState
                
                % to split, we increase the threshold for the duration of this
                % state and re-run the classification at the same time-scale.
                % Re-running instead of moving on avoids missing better peaks
                % to the right or left that were hidden inside the initial mask
                threshMap(stateMask1D) = min(1, threshMap(stateMask1D) + threshIncr);
                
                if all(threshMap(stateMask1D) == 1)
                    % this means we've exhausted the thresholding and the
                    % state is still not stable enough.
                    threshMap(stateMask1D) = Inf;
                else
                    incr = false;
                end
            else
                
                % take this state off the search grid:
                threshMap(stateMask1D) = Inf;
                
                kState = kState + 1;
                state(iTr).tm(kState, 1) = tmLo(max(1, dt(1))) + pad;
                state(iTr).tm(kState, 2) = tmLo(min(numel(tmLo), dt(2))) - pad;
                
            end
        end
        
        if incr
            pk = pk - 1; % move down one time-scale
        end
    end
   
    % STATE POSTPROCESSING:
     
    % now we go through every state and check whether
    % significance of the fit was lost for a longer period (Note that we
    % took the sig already into account above by assigning zero weights for
    % non significant fits when creating the cluster matrix, but
    % states could still be classified around gaps if the circular
    % stability was sufficiently high.)
    
    killState = false(1, size(state(iTr).tm,1));
    for iState = 1:size(state(iTr).tm,1)
        
        stateTimeWindow = (tmLo >= state(iTr).tm(iState,1) & tmLo <= state(iTr).tm(iState,2));
        sigMask = ~sigLo(:,iTr)' & stateTimeWindow;
        
        lbl = bwlabel(sigMask);
        killGaps = find(arrayfun(@(x) sum(lbl == x)>interpDur, 1:max(lbl)));
        
        if isempty(killGaps)
            continue
        end
        
        killGaps = arrayfun(@(x) [find(lbl == x,1,'first') find(lbl == x,1,'last')], killGaps, 'UniformOutput', false);
        killGaps = reshape([find(stateTimeWindow,1,'first') cell2mat(killGaps) find(stateTimeWindow,1,'last')], 2,[]);
        % -> now a 2xNumGaps matrix with on- and offsets of new states in
        % the first dim
        
        state(iTr).tm = cat(1, state(iTr).tm, tmLo(killGaps)'); % add new states to the end of the list
        killState(iState) = true; % mark this state for execution
        
        
    end
    
    state(iTr).tm(find(killState),:) = []; % delete previous longer versions
    
    
    % Now we sort all the states by their onset:
    [~, idx] = sort(state(iTr).tm(:,1));
    state(iTr).tm = state(iTr).tm(idx,:);

    
    
    % To adjust for temporal inaccuracies in our classification
    % (which is based on crude thresholding), we look for rapid transitions
    % in close proximity to state on- or offset, and adjust accordingly.
    isOK = false;
    updateShift = true;
    
    while ~isOK
        
        for iState = 1:size(state(iTr).tm,1)
            for iSide = 1:2 % onset/offset
                s2tDist = abs(rTrans{iTr} - state(iTr).tm(iState,iSide)); % temporal distance to each transition
                if ~any(s2tDist <= shift2TransMaxTime)
                    % none close enough
                    continue
                end
                tagme = find(s2tDist == min(s2tDist(s2tDist <= shift2TransMaxTime))); % pick the one closest
                
                if numel(tagme)==2
                    % if two have the exact same distance, we pick
                    % the first if we're at an offset, and the second if
                    % we're at an onset
                   tagme = tagme(1+(iSide==1)); 
                end
                    
                state(iTr).tm(iState,iSide) = rTrans{iTr}(tagme) - ((iSide*2)-3)*state2TransPad;
            end
        end
        
        updateShift = false;
        
        % Next, we look for overlapping states, delete the contended
        % period from both states and repeat the shifting to the nearest
        % transition if applicable
        for iState = 1:size(state(iTr).tm,1)
            olIdx = iState + find(state(iTr).tm(iState+1:end,1) < state(iTr).tm(iState,2), 1);
            if ~isempty(olIdx)
                newOff = state(iTr).tm(olIdx,1) - state2TransPad;
                newOn = state(iTr).tm(iState,2) + state2TransPad;
                state(iTr).tm(iState,2) = newOff;
                state(iTr).tm(olIdx,1) = newOn;
                updateShift = true;
            end
        end
        
        isOK = ~updateShift; % exit if no more update
        
    end
    
    
    % sort out states that are too short after all processing steps:
    state(iTr).tm(diff(state(iTr).tm,[],2) < p.Results.StateMinDur*1e-3, :) = [];
    
%     % sort out states that are not based on enough significant individual fits:
%     killState = false(1, size(state(iTr).tm,1));
%     for iState = 1:size(state(iTr).tm,1)
%         stateTimeWindow = (tmLo >= state(iTr).tm(iState,1) & tmLo <= state(iTr).tm(iState,2));
%         killState(iState) =  mean(~wav.sig(stateTimeWindow,iTr)) > sigFitThresh;
%     end
%     state(iTr).tm(killState,:) = [];
    
    
    % finally, we extract parameters for each state:
    wdFull = exp(1i.*wav.wavDir(:,iTr)); % retrieve the full wave dir, at original sampling rate & not masked by significance as above
    for iState = 1:size(state(iTr).tm,1)
        stateTimeWindow = (wav.t >= state(iTr).tm(iState,1) & wav.t <= state(iTr).tm(iState,2));
        rValid = mean(wdFull(stateTimeWindow & wav.sig(:,iTr)'));
        rInvalid = mean(wdFull(stateTimeWindow & ~wav.sig(:,iTr)'));
        state(iTr).dir(iState) = angle(rValid); % dir of mean angle (valid only)
        state(iTr).mag(iState) = abs(rValid); % magnitude of mean angle (valid only)
        state(iTr).dir_invalid(iState) = angle(rInvalid); % dir of mean angle (invalid only)
        state(iTr).mag_invalid(iState) = abs(rInvalid); % magnitude of mean angle (invalid only)
        state(iTr).pvalid(iState) = mean(wav.sig(stateTimeWindow,iTr)); 
    end
    

    %%
    
end

delete(wb);




%%
% output parameters:
param = p.Results;



end
