%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

% Computes pairs of nodes with a certain distance
%
% Input: 
%     neighborhoodSize: Defines the distance between node pairs in neighborhoodUnit
%     neighborhoodUnit: Either 'edgeCount' or 'meanEdgeLength'
%
% Output:
%     pairs: Nx2 node indices

function pairs = computeNodePairs(mesh, neighborhoodSize, neighborhoodUnit)

TR = triangulation(mesh.e, mesh.p);
edg = TR.edges;

switch neighborhoodUnit
    case 'edgeCount'
        dist = distances(graph(edg(:,1), edg(:,2)));
        d = 1;
    case 'meanEdgeLength'
        dist = computeSurfaceGeodesicMatrix(mesh);
        d = mean(sqrt(sum((mesh.p(edg(:,1),:)-mesh.p(edg(:,2),:)).^2,2)));
    otherwise
        error('Unknown neighborhoodUnit ''%s''.', neighborhoodUnit);
end

pairs = NaN(20*mesh.nop,2);
numPairs = 0;
for i = 1:mesh.nop
    neighs = find(dist(:,i) > (neighborhoodSize-0.5)*d & dist(:,i) < (neighborhoodSize+0.5)*d);
    numNeighs = numel(neighs);
    pairs(numPairs+1:numPairs+numNeighs,:) = [repmat(i,numNeighs,1) neighs];
    numPairs = numPairs + numNeighs;
end
pairs(isnan(pairs(:,1)),:) = [];
pairs = unique(sort(pairs,2), 'rows');

end