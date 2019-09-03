%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function [actTimesMat,Gx,Gy,Gz] = globalActTimes_prepare(mesh, pairs, weights)

if nargin < 3
    weights = ones(size(pairs,1),1);
end

% This is equivalent to the commented out code with for-loops below:
i = [repmat((1:size(pairs,1))',2,1); repmat(size(pairs,1)+1,mesh.nop,1)];
j = [pairs(:,1); pairs(:,2); (1:mesh.nop)'];
v = [ones(size(pairs,1),1); -ones(size(pairs,1),1); ones(mesh.nop,1)];
M = sparse(i, j, v, size(pairs,1)+1, mesh.nop);
M = [weights; 1] .* M;

% D = zeros(size(pairs,1), mesh.nop);
% for i = 1:size(pairs,1)
%     D(i,pairs(i,1)) = 1;
%     D(i,pairs(i,2)) = -1;
% end
% I = ones(1,mesh.nop);
% M = sparse(cat(1,D,I));

W = spdiags([weights(:); 1], 0, size(pairs,1)+1, size(pairs,1)+1);

actTimesMat = (M'*W*M)\M'*W;

if nargout > 1
    [Gx,Gy,Gz] = Gradient(mesh);
end

end