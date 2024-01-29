function [alphaP, alphaC] = generateAlphaSinusoids(fp, fc, sinrange, params)
% generates two sinusoidal luminance modulation timeseries
% one for periphery, one for center.
% fc, fp = frequencies in hertz
% sinrange: between 0 and 1, what should the luminance range be.
% sinrange of 1 = luminance ranges between 0 and 2*perceived
% mean luminance is always 1. That corresponds to the target perceived luminance

domain = 0:params.real_IFI:params.trial_length_max-params.real_IFI; % total time in s, at increments of 1 frame
fpRad = fp*2*pi; % periphery
fcRad = fc*2*pi; % center
% luminance modulation
alphaP = sinrange*sin(fpRad.*domain) + 1;
alphaC = sinrange*sin(fcRad.*domain) + 1;
% alpha transparency modulation
%alphaP = 0.5*sinrange*sin(fpRad.*domain) + 0.5;
%alphaC = 0.5*sinrange*sin(fcRad.*domain) + 0.5;
end