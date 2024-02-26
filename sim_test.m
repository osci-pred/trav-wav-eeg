clearvars

% simulate source signal:
simArgs = {...
    'Dur', 3,...
    'Onset', 1,...
    };

% simArgs = {...
%     'Dur', 6,...
%     'Frequency', 10,...
%     'On', [1 3],...
%     };

nTrials = 500;

for iTrial = 1:nTrials
    [src(:,1,iTrial), t] = sim_getSource('ERP', simArgs{:});
end


% project to 3D EEG elec positions:
projArgs = {...
    'SourceAlignElec', 'POz',...
    'NumNoiseSources', 8,...
    'snrRange', [1.2 3],...
    'Plot', true,...
    };

[eeg, t, elecLbl] = sim_eegProjection(src, t, projArgs{:});



% do the wave fit:
fitParam = {...
    'Frequency', [7 13],...
    'WindowSize', 100,... % in ms
    'NumStepsSpatFreq', 30,...
    'NumStepsWaveDir', 60,...
    };

wav = tw_fitPlaneEEG(eeg, t, elecLbl, fitParam{:});


isfw = @(x) abs(x+pi/2) < 0.5;
isbw = @(x) abs(x-pi/2) < 0.5;

%% plot 
figure
tiledlayout(3,1)

nexttile(1)
plot(t, squeeze(src));

nexttile(2)
plot(t, squeeze(eeg(:,:,1)));

nexttile(3)
plot(wav.t, mean(isfw(wav.wavDir),2), 'r');
hold on
plot(wav.t, mean(isbw(wav.wavDir),2), 'b');
legend({'FW' 'BW'})
ylabel('Prob.')
xlabel('Time [sec]')

