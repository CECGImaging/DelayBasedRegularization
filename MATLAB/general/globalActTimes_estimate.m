%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function actTimes = globalActTimes_estimate(X, pairs, actTimesMat, gaussSigma, mode, Gx, Gy, Gz)

%% Upsample and filter original signal

upsampFactor = 10;
sig = interp1(0:size(X,2)-1, X', 0:1/upsampFactor:size(X,2)-1)';

if gaussSigma > 0
    gaussSigma = upsampFactor * gaussSigma;
    bt = sqrt(log(2))/(2*pi*gaussSigma);
    span = round(6*gaussSigma/2)*2;
    gaussCoeffs = gaussdesign(bt, span, 1);
    sig = padarray(sig', span, 'replicate');
    sig = filtfilt(gaussCoeffs, 1, sig);
    sig = sig(1+span:end-span,:)';
end

%% Define derivative signal for cross-correlation

switch mode
    case 'spatial'
        derivSig = sqrt((Gx*sig).^2+(Gy*sig).^2+(Gz*sig).^2);
    case 'temporal'
        derivSig = (sig(:,3:end)-sig(:,1:end-2))/2;
    case 'spatiotemporal'
        gradSig = sqrt((Gx*sig).^2+(Gy*sig).^2+(Gz*sig).^2);
        diffSig = (sig(:,3:end)-sig(:,1:end-2))/2;
        derivSig = gradSig(:,2:end-1) .* diffSig;
    otherwise
        error('Unknown mode ''%s''.', mode);
end

%% Compute delays

delays = NaN(size(pairs,1),1);
parfor i = 1:size(pairs,1)
    [xc,lag] = xcorr(derivSig(pairs(i,1),:), derivSig(pairs(i,2),:));
    [~,ind] = max(xc);
    delays(i) = lag(ind);
end
delays = [delays; 0];

%% Compute actTimes using least-squares regression

actTimes = actTimesMat * delays;

%% Determine constant offset of activation times from mean of aligned signals

derivSigAligned = delayseq(derivSig', -actTimes);
[~,offset] = max(mean(derivSigAligned,2));
actTimes = (actTimes+offset)/upsampFactor;

end