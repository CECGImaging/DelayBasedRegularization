%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function [Gx,Gy,Gz] = Gradient(mesh, AtNodes)
% Calculates operator matrices for the surface gradient 
% of a triangular mesh using linear shape functions.
% Written by Steffen Schuler in March 2017,
% based on libigl's grad_tri() function.
% Input:  Mesh in the BEM-library format.
%         Set AtNodes to 0, if you want to compute the gradient at the
%         triangle centroids instead of at the nodes. Default is 1.
% Output: Gradient operator matrices for x,y,z components of the gradient.

if nargin < 2
	AtNodes = 1;
end

nop = mesh.nop;
noe = mesh.noe;

eperp21 = zeros(noe,3);
eperp13 = zeros(noe,3);

for i = 1:noe
    i1 = mesh.e(i,1);
    i2 = mesh.e(i,2);
    i3 = mesh.e(i,3);
    
    % compute triangle edge vectors
    v32 = mesh.p(i3,:) - mesh.p(i2,:);
    v13 = mesh.p(i1,:) - mesh.p(i3,:);
    v21 = mesh.p(i2,:) - mesh.p(i1,:);
    
    % compute normal vector
    n = cross(v32, v13);
    dblA = norm(n);
    n = n / dblA;
    
    % rotate each edge vector 90 degrees around normal 
    eperp21(i,:) = cross(n, v21);
    eperp21(i,:) = eperp21(i,:) / norm(eperp21(i,:));
    eperp21(i,:) = eperp21(i,:) * norm(v21) / dblA;
    eperp13(i,:) = cross(n, v13);
    eperp13(i,:) = eperp13(i,:) / norm(eperp13(i,:));
    eperp13(i,:) = eperp13(i,:) * norm(v13) / dblA;
end

% row indices
rs = zeros(1,3*4*noe);
k = 1;
for r = 1:3
    for j = 1:4
        for i = ((r-1)*noe+1):r*noe 
            rs(k) = i;
            k = k+1;
        end
    end
end

% column indices
cs = zeros(1,3*4*noe);
k = 1;
for r = 1:3
    for i = 1:noe; cs(k) = mesh.e(i,2); k = k+1; end
    for i = 1:noe; cs(k) = mesh.e(i,1); k = k+1; end
    for i = 1:noe; cs(k) = mesh.e(i,3); k = k+1; end
    for i = 1:noe; cs(k) = mesh.e(i,1); k = k+1; end
end

% values
vs = zeros(1,3*4*noe);
k = 1;
for i = 1:noe; vs(k) =  eperp13(i,1); k = k+1; end
for i = 1:noe; vs(k) = -eperp13(i,1); k = k+1; end
for i = 1:noe; vs(k) =  eperp21(i,1); k = k+1; end
for i = 1:noe; vs(k) = -eperp21(i,1); k = k+1; end
for i = 1:noe; vs(k) =  eperp13(i,2); k = k+1; end
for i = 1:noe; vs(k) = -eperp13(i,2); k = k+1; end
for i = 1:noe; vs(k) =  eperp21(i,2); k = k+1; end
for i = 1:noe; vs(k) = -eperp21(i,2); k = k+1; end
for i = 1:noe; vs(k) =  eperp13(i,3); k = k+1; end
for i = 1:noe; vs(k) = -eperp13(i,3); k = k+1; end
for i = 1:noe; vs(k) =  eperp21(i,3); k = k+1; end
for i = 1:noe; vs(k) = -eperp21(i,3); k = k+1; end

% create sparse gradient operator matrix
G = sparse(rs,cs,vs);
Gx = G(1:noe,:);
Gy = G(noe+1:2*noe,:);
Gz = G(2*noe+1:end,:);

if AtNodes
    % create sparse matrix T transforming gradient at triangle centroids
    % into grad at nodes by averaging gradients of neighboring triangles
    rows = zeros(1,10*nop);
    cols = zeros(1,10*nop);
    vals = zeros(1,10*nop);
    k = 1;
    for i = 1:nop
        notri = mesh.ntri_n(i);
        tri = mesh.ntri(i,1:notri);

        range = k:k+notri-1;
        rows(range) = repmat(i, 1,notri);
        cols(range) = tri;
        vals(range) = repmat(1/notri, 1,notri);

        k = k+notri;
    end
    ind = find(cols==0);
    cols(ind) = [];
    rows(ind) = [];
    vals(ind) = [];
    
    T = sparse(rows,cols,vals);

    Gx = T * Gx;
    Gy = T * Gy;
    Gz = T * Gz;
    % spy(Gx);
end

end
