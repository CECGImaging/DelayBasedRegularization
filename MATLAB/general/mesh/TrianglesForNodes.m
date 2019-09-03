function [tr,sides,n]=TrianglesForNodes(triangles,non)
% function [tr,sides,n]=TrianglesForNodes(triangles,non)
% Finds, which triangles belong to the neighborhood of each node
% tr = list of triangles for each node
% sides = index of the node in the triangle (1, 2, or 3)
% n = number of triangles for each node

maxt=20;%maximum number of triangles for one node
tr=zeros(non,maxt);
sides=zeros(non,maxt);
n=zeros(non,1);
not=size(triangles,1);
% count=zeros(non,1);
e1=triangles(:,1);
e2=triangles(:,2);
e3=triangles(:,3);

for I=1:not
    n(e1(I))=n(e1(I))+1;
    tr(e1(I),n(e1(I)))=I;
    sides(e1(I),n(e1(I)))=1;
    
    n(e2(I))=n(e2(I))+1;
    tr(e2(I),n(e2(I)))=I;
    sides(e2(I),n(e2(I)))=2;
    
    n(e3(I))=n(e3(I))+1;
    tr(e3(I),n(e3(I)))=I;
    sides(e3(I),n(e3(I)))=3;

end
nmax=max(n);
tr=tr(:,1:nmax);
sides=sides(:,1:nmax);
