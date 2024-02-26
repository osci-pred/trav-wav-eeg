function [fwMax,bwMax,midline,freq]=tw_wavesHunterAllFreqs(data,samplingRate)

    %data is a matrix with [electrode, timePoints] as dimension. The
    %electrode dimension should be odd (it simplifies the identificaton of
    %the midline)
    % it returns three spectra for fw, bw waves, and the midline (standing
    % waves)

    durationSignal=size(data,2)/samplingRate; %it's in second.
    dF = 1/durationSignal;                   
    fx = -samplingRate/2:dF:samplingRate/2-dF;         
    fy=(size(data,1)/2)*linspace(-1,1,size(data,1));
    sbj2DFFT=abs(fftshift(fft2(data)));
    
    fwMax=max(sbj2DFFT(fy<0, fx>=1 & fx<=45));
    bwMax=max(sbj2DFFT(fy>0, fx>=1 & fx<=45));
    midline=sbj2DFFT(fy==0, fx>=1 & fx<=45);
    
    freq = fx(fx>=1 & fx<=45);
end

    