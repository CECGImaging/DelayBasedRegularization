%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function L = LaplaceBeltrami(mesh)
% Calculates the operator matrix for the surface Laplacian
% of a triangular mesh (Laplace-Beltrami operator).
% Written by Steffen Schuler in March 2017,
% based on libigl.
% Input:  Mesh in the BEM-library format.
% Output: Laplace operator matrix.

[Gx,Gy,Gz] = Gradient(mesh, 0);
G = [Gx; Gy; Gz];
T = spdiags(repmat(mesh.a,3,1), 0, 3*mesh.noe, 3*mesh.noe);
L = -G' * T * G;

end
