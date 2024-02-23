function [src, t] = sim_getSource(type, varargin)
% get dipole sources for simulation

switch type
    case 'Sine'
        [src, t] = sineWave(varargin{:});
    case 'ERP'
        [src, t] = erp(varargin{:});
end

end
%%
function [out, t] = sineWave(varargin)
p = inputParser;
p.addParameter('numSources', 1); % number of sources
p.addParameter('Dur', 1); % duration in sec
p.addParameter('Frequency', 10); 
p.addParameter('PhaseOffset', 0); % in radians
p.addParameter('Baseline', 0); % corresponding to off state
p.addParameter('Offset', 0); % corresponding to mean of sine
p.addParameter('On', [-Inf Inf]); % on times (ouside is baseline)
p.addParameter('Scale', 20); % Peak to trough
p.addParameter('Fs', 160); % sampling frequency Hz

if strcmpi(varargin{1}, 'help')
    p.parse();
    t = [];
    out = p.Results;
    return
else
    p.parse(varargin{:});
end

on = p.Results.On;
if ~iscell(on)
    on = {on};
end

val = p.Results.Scale;

t = 1/p.Results.Fs:1/p.Results.Fs:p.Results.Dur;

out = ones(1,numel(t)).*p.Results.Baseline;
for ievent = 1:numel(on)
    out(rngfn(t, on{ievent})) = 1;
end

out = out .* (sin((2*pi*p.Results.Frequency).*t + p.Results.PhaseOffset) .* val/2 + val/2 + p.Results.Offset);
out = out';

end
%%
function [out, t] = erp(varargin)
p = inputParser;
p.addParameter('Dur', 2, @(x) x > 1); % duration sec
p.addParameter('Baseline', 0); % corresponding to off state
p.addParameter('Onset', 0.2); % on times (ouside is baseline)
p.addParameter('Scale', 20); % Peak
p.addParameter('Fs', 160); % sampling frequency Hz

if strcmpi(varargin{1}, 'help')
    p.parse();
    t = [];
    out = p.Results;
    return
else
    p.parse(varargin{:});
end

t = 1/p.Results.Fs:1/p.Results.Fs:p.Results.Dur;
t = t - p.Results.Onset;

% ERP template is fixed, we just scale to size:
modHi = normpdf(t, 0.15, 0.05);
sinHi = sin(2*pi*(5.6).*t + 0.8);
modLo = 1.1.*normpdf(t, 0.35, 0.15);
sinLo = sin(2*pi*2.9.*t);

template = modHi.*sinHi + modLo.*sinLo;
template = template ./ max(abs(template));

out = template .* p.Results.Scale + p.Results.Baseline;
out = out';

end
%%
function out = rngfn(vec, rng)
% check where array is in range
out = vec >= rng(1) & vec < rng(2);
end
