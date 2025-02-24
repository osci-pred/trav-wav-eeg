clearvars -except dataA3
isubj = 13;

data = dataA3(isubj);

sel = data.dynamicEEG;

roi = tw_getROI(data.channelInfo,...
    {'I' 'O' 'PO' 'P' 'CP' 'C' 'FC'},... % Inion to Fronto Central
    0:2); % z,1,2


wav = tw_fitPlaneEEG(sel, data.timepoints, data.channelInfo, 'Frequency', [7 13], 'WindowSize', 50, 'ROI', roi);

state = tw_stateSwitchAnalysis(wav);


%% plot a single trial
figure
iTr = 50;
hold on
subplot(211)
midLine = tw_getMidLine(data.channelInfo, 1:numel(data.channelInfo));
contourf(wav.t, 1:numel(midLine), sel(midLine, :,iTr), 64, 'EdgeColor', 'none')
xlim([min(wav.t) max(wav.t)])
yticks([1 numel(midLine)])
yticklabels(data.channelInfo(midLine([1 end])))
colormap(brewermap(64, '*greens'))

subplot(212)
plot(wav.t, wav.wavDir(:,iTr));
stlbl = {'fw' 'bw'};
cols = {'b' 'r'}; % for patches
for ist = 1:2
    st = stlbl{ist};
    mask = state.(st){iTr};
    
    for i = 1:size(mask,1)
        patch(mask(i, [1 2 2 1]), [-1 -1 1 1].*pi, [1 1 1 1], cols{ist}, 'EdgeColor', 'none');
        alpha 0.3
    end
end
xlim([min(wav.t) max(wav.t)])

%% plot histogram of (stable) state durations
figure

binEdge = state.param.StateMinDur*1e-3:0.05:1.5;
statedur.bins = binEdge(1:end-1) + 0.5*diff(binEdge(1:2));

for st = {'fw' 'bw'}
statedur.(char(st)) = histcounts(diff(cat(1,state.(char(st)){:}), [],2), binEdge, 'Normalization', 'probability');
end

hold on
plot(statedur.bins, statedur.fw);
plot(statedur.bins, statedur.bw);
xline(state.param.StateMinDur*1e-3, 'k:')
xlim([0 1.5])
xlabel('Duration [secs]')
ylabel('Prob.')
title('Duration of stable states')
%% plot probability of stable FW and BW state over time:

figure
subplot(211)
hold on
plot(wav.t, wav.pFitFW, 'b')
plot(wav.t, wav.pFitBW, 'r')
legend({'FW' 'BW'}, 'AutoUpdate', 'off')
title('Original classification (significance only)')
ylim([-0.1 1.1])
yline(1, 'k:')
yline(0, 'k:')
ylabel('P (state)')

subplot(212)
hold on
plot(wav.t, 1-nanmean(state.maskedBoolFW, 2), 'r')
ylabel('P (stable state is BW)')
ylim([-0.1 1.1])
yyaxis right
plot(wav.t, mean(~isnan(state.maskedBoolFW),2), 'k');
ylabel('P (any stable state)')
ylim([-0.1 1.1])

yline(1, 'k:')
yline(0, 'k:')

yline(0.5, 'k:')

xlabel('Time [secs]')

ax = gca;
ax.YAxis(2).Color = 'k';