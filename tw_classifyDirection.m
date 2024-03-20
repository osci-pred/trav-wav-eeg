function [fw, bw] = tw_classifyDirection(wavDir, winSize)
fw = abs(wavDir+pi/2) < winSize;
bw = abs(wavDir-pi/2) < winSize;
end