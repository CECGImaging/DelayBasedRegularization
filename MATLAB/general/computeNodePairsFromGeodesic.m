%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

% Computes pairs of nodes and corresponding distance weigths from a
% geodesic distance matrix
%
% Input: 
%     geodesic:    numNodes x numNodes geodesic distance matrix
%     maxDistance: maximum distance between pairs of nodes
%                  (same unit as geodesic distances)
%
% Output:
%     pairs:   Nx2 node indices
%     weights: Nx1 cosine distance weights. The weight is 1 at distance 0,
%              0.5 at maxDistance/2 and 0 at maxDistance or larger

function [pairs,weights] = computeNodePairsFromGeodesic(geodesic, maxDistance)

numNodes = size(geodesic,1);
pairs = NaN(100*numNodes,2);
numPairs = 0;
for i = 1:numNodes
    neighs = find(geodesic(:,i) < maxDistance);
    numNeighs = numel(neighs);
    rowInd = numPairs+1:numPairs+numNeighs;
    pairs(rowInd,:) = [repmat(i,numNeighs,1) neighs];
    numPairs = numPairs + numNeighs;
end
pairs(numPairs+1:end,:) = [];
pairs(diff(pairs,1,2)==0,:) = [];
pairs = unique(sort(pairs,2), 'rows');
numPairs = size(pairs,1);

distanceWeight = @(x) 0.5*(cos(pi/maxDistance*x)+1);
weights = NaN(numPairs,1);
for i = 1:numPairs
    weights(i) = distanceWeight(geodesic(pairs(i,1),pairs(i,2)));
end

end