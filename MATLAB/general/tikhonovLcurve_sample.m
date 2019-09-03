%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function [logRho,logEta,lambda,curvature] = tikhonovLcurve_sample(B, A, R, zeroMeanNullspaceRemoval, limits, numSamples, smoothingParam)

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

AA = A'*A;
RR = R'*R;
AB = A'*B;

if zeroMeanNullspaceRemoval
    AA = AA + repmat(mean(abs(AB(:)))/size(A,2), size(A,2));
end

lambda = logspace(limits(1), limits(2), numSamples);
solNorm = NaN(numSamples,1); % solution norm
resNorm = solNorm;           % residual norm

for i = 1:numel(lambda)
    X = (AA + lambda(i)*RR) \ AB;
    resNorm(i) = norm(A*X-B, 'fro');
    solNorm(i) = norm(R*X, 'fro');
end

logRhoOrig = log10(resNorm);
logEtaOrig = log10(solNorm);
ft = fittype('smoothingspline');
opts = fitoptions('Method', 'SmoothingSpline');
opts.SmoothingParam = smoothingParam;
fitresult = fit(logRhoOrig, logEtaOrig, ft, opts);

logRho = linspace(logRhoOrig(1), logRhoOrig(end), 10*numel(logRhoOrig));
logEta = feval(fitresult, logRho);
lambda = spline(logRhoOrig, lambda, logRho);

[d1,d2] = differentiate(fitresult, logRho);
curvature = d2./(1+d1.^2).^(3/2);

end
