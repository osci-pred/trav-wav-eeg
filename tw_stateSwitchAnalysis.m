function out = tw_stateSwitchAnalysis(wav, varargin)
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
% 'StateMinDur'             min duration of state (in ms).
% 'StateInterpMaxDur'       max duration of gap (non-significant fit or leaving dir window) to
%                           interpolate over during stable state.
% 'StateTransTagTol'        tolerance for the tmp. distance (in ms) between a transition
%                           and state onset/offset to be tagged as belonging together

p = inputParser();
p.addParameter('TransThresh', 0.2);
p.addParameter('TransMaxDur', 100);
p.addParameter('TransMinDist', pi/4);
p.addParameter('StateMinDur', 100);
p.addParameter('StateInterpMaxDur', 100);
p.addParameter('StateTransTagTol', 50);

p.parse(varargin{:});

%% Step 1: detect rapid transitions in the wave direction

tStep = diff(wav.t(1:2));
sr = 1/tStep;
nTr = size(wav.wavDir,2);

maxVelEvDur = p.Results.TransMaxDur*1e-3 * sr; % max high-vel event dur in samples

vel = nan(size(wav.wavDir));
out.trans = cell(1,nTr);

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
    out.trans{iTr} = arrayfun(@(evId) min(wav.t(abs(vel(:,iTr))'.*(ev == evId) == max(abs(vel(:,iTr))' .* (ev == evId)))), find(inclEv));
    out.trans{iTr} = out.trans{iTr} - 0.5*tStep; % account for left-side diff above
%     
%         % plot for sanity check:
%         plot(wav.t, wav.wavDir(:,iTr));
%         hold on
%         for i = find(inclEv)
%             tst = wav.t(find(ev == i, 1, 'first')-1);
%             tend = wav.t(find(ev == i,1, 'last')+1);
%             patch([tst tend tend tst], [-1 -1 1 1].*pi, [1 1 1 1], 'EdgeColor', 'none');
%             alpha 0.3
%         end
%         xline(out.trans, 'k:')
%         keyboard
end


%% Step 2: detect and interpolate stable states of FW and BW waves


% for now, we use the classification done in tw_fitPlaneEEG for this, but
% maybe this can be improved (re-classify based on peaks for each subject?)

% wav.fw and wav.bw are boolean indicating classified states
% wav.sig is boolean indicating whether the fit was considered significant

stateDirs = {'fw' 'bw'};
minDur = p.Results.StateMinDur*1e-3 * sr; %min dur in samples
interpDur = p.Results.StateInterpMaxDur*1e-3 * sr; % in samples

transCoincidePad = 0.05; % in secs, how far to each side we check whether a state-break coincides with a transition event

for ist = 1:2
    
    st = stateDirs{ist};
    out.(st) = cell(1,nTr);
    
    for iTr = 1:nTr
        
        
        stateVec = wav.(st)(:,iTr) & wav.sig(:,iTr);
        lbl = bwlabel(~stateVec); % label OFF-states to identify state-breaks
        evDurCrit = arrayfun(@(id) sum(lbl==id) <= interpDur, unique(lbl(lbl>0))); % max dur
        
        hasTransCrit = arrayfun(@(id) any(any(abs(wav.t(lbl==id)-out.trans{iTr}') < transCoincidePad)), unique(lbl(lbl>0)));
        % -> checks whether event coincides with a transition event
        
        % interpolate all the events satisfying both criteria:
        stateVec(ismember(lbl, find(~hasTransCrit & evDurCrit))) = true;
        
        % relabel to get ON-states:
        lbl = bwlabel(stateVec);
        
        % apply min duration criterion:
        inclEv = arrayfun(@(evID) sum(lbl == evID) >= minDur, unique(lbl(lbl>0)));
        stateVec(~ismember(lbl, find(inclEv))) = 0;
        lbl = bwlabel(stateVec); % get final state labels
        
        
        out.(st){iTr} = cell2mat(arrayfun(@(evID)...
            wav.t([find(lbl==evID, 1, 'first'), find(lbl==evID, 1, 'last')]),...
            unique(lbl(lbl>0)), 'UniformOutput', false));
    end
    
end


% % plot for sanity check:
% plot(wav.t, wav.wavDir(:,iTr));
% hold on
% for i = unique(lbl(lbl>0))'
%     tst = wav.t(find(lbl == i, 1, 'first'));
%     tend = wav.t(find(lbl == i,1, 'last'));
%     patch([tst tend tend tst], [-1 -1 1 1].*pi, [1 1 1 1], 'EdgeColor', 'none');
%     alpha 0.3
% end
% xline(out.trans{iTr}, 'k:')



%% add information about neighbouring states to the transition times:

sides = {'pre' 'post'};

for iTr = 1:nTr
    out.tag{iTr} = zeros(numel(out.trans{iTr}),2);
    
    for iSide = 1:2
        
        for iTrans = 1:numel(out.trans{iTr})
            
            for ist = 1:2
                
                st = stateDirs{ist};
                
                if isempty(out.(st){iTr})
                    continue
                end
                
                distCheck = abs(out.trans{iTr}(iTrans)-out.(st){iTr}(:,iSide)) <= p.Results.StateTransTagTol*1e-3;
                
                if any(distCheck)
                    if out.tag{iTr}(iTrans,iSide) ~= 0
                        error('End of the world') % this shouldn't happen
                    end
                    out.tag{iTr}(iTrans,iSide) = ist; % save index [1/2] for {'fw' 'bw'}
                end
            end
            
        end
    end
end
out.tagLabels = {'fw' 'bw'};


%% return continuous state array for each trial, as in wav struct before, 
% but this time masked by nans wherever state is not one of FW or BW:

rangefn = @(x, lms) x >= lms(1) & x <= lms(2);
for ist = 1:2
    st = stateDirs{ist};
    tmp.(st) = zeros(size(wav.fw));
    for iTr = 1:nTr
        
        for iEv = 1:size(out.(st){iTr},1)
            tmp.(st)(rangefn(wav.t, out.(st){iTr}(iEv,:)),iTr) = 1;
        end
    end
end

idx = ~(tmp.fw | tmp.bw);
out.maskedBoolFW = tmp.fw;
out.maskedBoolFW(idx) = nan;
% maskedBoolFW now contains entries (0: BW, 1: FW) wherever either state
% was classified, nans otherwise

% save parameters in output:
out.param = p.Results;



end
