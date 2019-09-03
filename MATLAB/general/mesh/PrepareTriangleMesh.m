function mesh=PrepareTriangleMesh(points,elements)
% function mesh=PrepareTriangleMesh(points, elements);
% Calculates normal vectors, areas, and other things needed in BEM transfer
% matrix generation
% mesh.p = nodes of the mesh
% mesh.e = element description; points to the nodes of the mesh
% mesh.nop = number of nodes
% mesh.noe = number of triangles
%
% mesh.n = triangle normals, length 2*TriangleArea
% mesh.un = triangle unit normals, 
% mesh.a = triangle areas
% mesh.mp = triangle midpoints
%
% The following ones are unnecessary for BEM matrix generation, but they
% are often of some use
% mesh.ntri_n = number of triangles belonging to each node
% mesh.ntri = triangles for each node
% mesh.ntri_s = indices of the nodes in each of the belonging triangles
% mesh.nn = unit normals of the nodes (approximation)

points=double(points);
elements=double(elements);

[mesh.n,mesh.un,mesh.a]=TriangleNormals(points,elements);
mesh.nop=size(points,1);
mesh.noe=size(elements,1);
[mesh.ntri,mesh.ntri_s,mesh.ntri_n]=TrianglesForNodes(elements,mesh.nop);
mesh.nn=NodeNormals(mesh.un,mesh.ntri,mesh.ntri_n);
mesh.mp=TriangleMidpoints(points,elements);
mesh.p=points;
mesh.e=elements;
