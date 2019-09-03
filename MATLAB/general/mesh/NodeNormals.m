function nnormals=NodeNormals(tnormals,ntri,ntri_n)
% function nnormals=NodeNormals(tnormals,ntri,ntri_n)
% Calculates unit normals for the nodes
% tnormals = normals of the triangles
% tnormals can be either unit normals (mesh.un) or normals with length 2A
% (mesh.n); in a nonregular mesh the results are somewhat different - see,
% which best suits your purpose.
non=size(ntri,1);
nnormals=zeros(non,3);
for I=1:non
    tempnormal=[0 0 0];
    for J=1:ntri_n(I)
        tempnormal = tempnormal + tnormals(ntri(I,J),:);
    end
    nnormals(I,:)=tempnormal/(sqrt(tempnormal*tempnormal'));
end
