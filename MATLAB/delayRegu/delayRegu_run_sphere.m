%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

addpath(genpath('../general'));

%% Config

inputDir = '../../ExampleData/sphere';
outputDir = '../../Reconstructions/sphere';

SNR = -10;                      % signal-to-noise ratio
numIter = 100;                  % number of iterations
numEdges = 2;                   % number of edges between node pairs
lcurveStartEndThresh = 0.15;    % threshold for defining the time interval used for the L-curve
lcurveBounds = [-5 0];          % bounds of logLambda for sampling the L-curve
actTimesSigma = 12;             % std for temporal Gaussian filter used in activation times estimation

% Nodes used for visualization of time courses
sensorNodes = [304 435 313 326 337 462 346 359 370]; % line along equator

%% General

load(sprintf('%s/geometry/heart_epi.mat', inputDir));
mesh = PrepareTriangleMesh(heart_epi.points, heart_epi.cells);

L = LaplaceBeltrami(mesh);

load(sprintf('%s/transferMat/transferTMV.mat', inputDir));
epiEndoProj = cat(1, speye(mesh.nop), speye(mesh.nop));
A = transferTMV * epiEndoProj;

load(sprintf('%s/signals/tmv.mat', inputDir));
tmv = double([repmat(tmv(:,1),1,50) tmv]);
tmv_epi = tmv(1:mesh.nop,:);

B = transferTMV * tmv;
rng(1);
B = addwhitenoise(B, SNR, [50 550]);

%% Determine node neighbor pairs and prepare activation times estimation

pairs = computeNodePairs(mesh, numEdges, 'edgeCount');

[actTimesMat,Gx,Gy,Gz] = globalActTimes_prepare(mesh, pairs);

%%
try

resultsDir = outputDir;
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end
figDir = sprintf('%s/fig', resultsDir);
if ~exist(figDir,'dir')
    mkdir(figDir);
end

diary(sprintf('%s/log.txt', resultsDir));

%% Estimate true activation times and visualize and store the truth

actTimes_true = globalActTimes_estimate(tmv_epi, pairs, actTimesMat, actTimesSigma, 'spatiotemporal', Gx, Gy, Gz);

fig = figure('Name','Truth', 'WindowStyle','docked');
delayRegu_visualize(actTimes_true, tmv_epi, mesh, sensorNodes);
saveas(fig, sprintf('%s/truth.fig', figDir));

results.TMV_true = tmv_epi;
results.AT_true  = actTimes_true;

%% Perform initial reconstruction using second order Tikhonov

% determine depolarization interval for L-curve from sum of BSP magnitudes
Bsum = movmean(sum(abs(B),1),50);
Bsum = (Bsum-min(Bsum))/(max(Bsum)-min(Bsum));
startInd = find(Bsum>lcurveStartEndThresh,1,'first');
endInd = find(Bsum(:,1+startInd:end)<lcurveStartEndThresh,1,'first');
endInd = endInd+startInd-1;

% plot normalized sum of BSP magnitudes
fig = figure;
set(fig, 'pos', get(fig,'pos').*[1 1 1.5 1]);
subplot(1,2,1)
plot(Bsum)
hold on
plot(startInd, Bsum(startInd), 'kx')
plot(endInd, Bsum(endInd), 'kx')
xlabel('Time')
title('Normalized sum of BSP magnitudes')

% determine initial lambdaL using L-curve
subplot(1,2,2)
lambdaL0 = tikhonovLcurve_corner(B(:,startInd:endInd), A, L, true, lcurveBounds, 50, 1, gca);
saveas(fig, sprintf('%s/Lcurve.fig', resultsDir));
close(fig);

X = tikhonov(B, A, L, lambdaL0, true);
residualNorm0 = norm(A*X-B,'fro');

X = X-mean(X(:));

% estimate activation times
actTimes = globalActTimes_estimate(X, pairs, actTimesMat, actTimesSigma, 'spatiotemporal', Gx, Gy, Gz);

% visualize results
fig = figure('Name','Reconstruction', 'WindowStyle','docked');
delayRegu_visualize(actTimes, X, mesh, sensorNodes);
saveas(fig, sprintf('%s/iter%03i.fig', figDir, 0));

% initialize and store results
results.TMV = NaN(size(A,2), size(B,2), numIter, 'single');
results.AT = NaN(size(A,2), numIter, 'single');
results.TMV(:,:,1) = X;
results.AT(:,1) = actTimes;

%% Find root of residual norm difference using secant method

% define factors used to decrease lambdaL throughout iterations
factorsLambdaL = linspace(1,0,numIter);

% history for extrapolation of logLambdaD, has to be initialized with NaNs
s = NaN(1,4);

tic
for i = 2:numIter
    
    iter = i-1;
    
    factorLambdaL = factorsLambdaL(i);
    fprintf('\niteration %i\t| factorLambdaL = %.2e\n', iter, factorLambdaL);
    
    lambdaL = factorLambdaL * lambdaL0;
    
    [X,s] = delayRegu_secant(s, residualNorm0, B, A, L, pairs, actTimes, lambdaL, X);
    X = X-mean(X(:));
    
    % estimate activation times
    actTimes = globalActTimes_estimate(X, pairs, actTimesMat, actTimesSigma, 'spatiotemporal', Gx, Gy, Gz);
    
    % visualize results
    figure(fig);
    delayRegu_visualize(actTimes, X, mesh, sensorNodes);
    saveas(fig, sprintf('%s/iter%03i.fig', figDir, iter));
    
    % store results
    results.TMV(:,:,i) = X;
    results.AT(:,i) = actTimes;
    
end
toc

save(sprintf('%s/results.mat', resultsDir), 'results', '-v7.3');

catch err
    disp(getReport(err, 'extended', 'hyperlinks', 'on'))
end

diary('off');
