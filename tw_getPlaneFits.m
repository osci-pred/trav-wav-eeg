function [phiPred, wvdir, a, b, xi] = tw_getPlaneFits(pos, numStepsWaveDir, numStepsSpatFreq, numCyclesMax)
% get predicted phases for all plane fits. NumStepsWaveDir / SpatFreq are
% the resolution in radial and circular dimensions, numCyclesMax is the
% maximum radius (i.e. max spatial freq.), as number of cycles across the
% covered electrode space

deltapos = sqrt((pos(:,1)-pos(:,1)').^2 + (pos(:,2)-pos(:,2)').^2);
d = max(deltapos(:)); % this is the max distance covered by the electrode space

xiMax = (numCyclesMax*2*pi)/d;

xiStep = xiMax/numStepsSpatFreq;
thetaStep = (2*pi)/numStepsWaveDir;

theta = -pi:thetaStep:pi-thetaStep; % circular parameter (-> wave direction)
a = (xiStep:xiStep:xiMax)' .* sin(theta);
b = (xiStep:xiStep:xiMax)' .* cos(theta);

a = a(:);
b = b(:);

xi = sqrt(a.^2 + b.^2); % radial parameter / spatial freq
wvdir = atan2(b,a); 
phiPred = a.*pos(:,1)' + b.*pos(:,2)'; % predicted phases
end