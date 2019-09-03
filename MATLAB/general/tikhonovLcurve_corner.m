%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function cornerLambda = tikhonovLcurve_corner(B, A, R, zeroMeanNullspaceRemoval, limits, numSamples, smoothingParam, ax)

if nargin < 7
    smoothingParam = 1-1e-6;
end
if nargin < 6
    numSamples = 50;
end
if nargin < 5
    limits = [-5 5];
end
if nargin < 4
    zeroMeanNullspaceRemoval = false;
end

[logRho,logEta,lambda,curvature] = tikhonovLcurve_sample(B, A, R, zeroMeanNullspaceRemoval, limits, numSamples, smoothingParam);

curvature(curvature<0) = 0;
curvature = curvature/max(curvature);
startInd = find(curvature>0.5, 1, 'first');
endInd = find(curvature(1+startInd:end)<0.5, 1, 'first');
endInd = endInd+startInd-1;

subLimits = log10([lambda(startInd) lambda(endInd)]);
[subLogRho,subLogEta,subLambda,subCurvature] = tikhonovLcurve_sample(B, A, R, zeroMeanNullspaceRemoval, subLimits, numSamples, smoothingParam);

[~,cornerInd] = max(subCurvature);
cornerLambda = subLambda(cornerInd);

rho = 10.^logRho;
eta = 10.^logEta;
subRho = 10.^subLogRho;
subEta = 10.^subLogEta;

if nargin == 8
    loglog(ax, rho, eta, 'k')
    hold on
    loglog(ax, subRho, subEta, 'r')
    loglog(ax, subRho(cornerInd), subEta(cornerInd), 'bx')
    title(ax, sprintf('L-curve lambda = %.2e', cornerLambda))
    xlabel(ax, 'Residual norm')
    ylabel(ax, 'Solution norm')
end

end
