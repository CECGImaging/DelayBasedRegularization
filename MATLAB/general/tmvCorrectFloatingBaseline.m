%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function [tmv,depolarCenter,repolarCenter] = tmvCorrectFloatingBaseline(tmv, pval)

if nargin < 2
    pval = 10;
end

% The spatial median abs. dev. is max, when the number of depolarized and 
% non-depolarized nodes is approx. the same. Before this depolarizarion 
% center, the baseline is defined using the lower percentile of TMVs in
% space, while the upper percentile is used behind.
% The second largest peak in the MAD corresponds to the center of
% repolarization, which is used to switch back from the upper to the lower
% percentile for defining the baseline.

% temporal centers using median abs. dev.
madSig = movmean(mad(tmv,1,1),30);
[~,peakLoc] = findpeaks(madSig, 'MinPeakProminence',0.1*max(madSig), 'MinPeakDistance', min(200,size(tmv,2)-2), 'sortstr','descend');

if numel(peakLoc) < 1
    % no depolarization center found
    % -> use lower percentile for all timesteps
    depolarCenter = size(tmv,2);
    repolarCenter = size(tmv,2);
elseif numel(peakLoc) < 2
    % only depolarization center found
    % -> switch once between lower and upper percentiles
    depolarCenter = peakLoc;
    repolarCenter = size(tmv,2);
else
    % de- and repolarization center found
    % -> switch twice between lower and upper percentiles
    peakLoc = sort(peakLoc(1:2));
    depolarCenter = peakLoc(1);
    repolarCenter = peakLoc(2);
end

lp = prctile(tmv,pval,1);       % lower percentile
up = prctile(tmv,100-pval,1);   % upper percentile

baseline = zeros(1,size(tmv,2));
baseline(1:depolarCenter-1) = lp(1:depolarCenter-1);
baseline(depolarCenter:repolarCenter-1) = lp(depolarCenter) - up(depolarCenter)+up(depolarCenter:repolarCenter-1);
baseline(repolarCenter:end) = lp(depolarCenter)-up(depolarCenter)+up(repolarCenter) - lp(repolarCenter)+lp(repolarCenter:end);

tmv = tmv - baseline;

end