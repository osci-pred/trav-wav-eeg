function out = sim_eegNoise(N, Fs)
out = rand(1,N);
[wt, f] = cwt(out, Fs);
wt = wt./abs(wt);
out = icwt(wt./sqrt(f));
end
