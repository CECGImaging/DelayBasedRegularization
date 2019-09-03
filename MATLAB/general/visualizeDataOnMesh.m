%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function visualizeDataOnMesh(mesh, data, pLimits, numSteps, viewAngles)

if nargin < 3
    pLimits = [0 100];
end
if nargin < 4
    numSteps = 20;
end
if nargin < 5
    viewAngles = [0 0];
end

trisurf(mesh.e, mesh.p(:,1), mesh.p(:,2), mesh.p(:,3), data,  'facecolor','interp', 'edgecolor','none');
axis equal;
camproj('persp');
set(gca,'visible','off');
caxis([prctile(data,pLimits(1)) prctile(data,pLimits(2))]);
colormap(jet(numSteps));
view(viewAngles(1), viewAngles(2));

end