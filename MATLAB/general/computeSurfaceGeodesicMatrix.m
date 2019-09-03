%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

% Relies on toolbox_graph from Gabriel Peyre:
% https://www.numerical-tours.com/installation_matlab/

function D = computeSurfaceGeodesicMatrix(mesh)

D = NaN(mesh.nop);
p = mesh.p;
e = mesh.e;
parfor i = 1:mesh.nop
    D(:,i) = perform_fast_marching_mesh(p, e, i);
end
D = (D+D')/2;

end