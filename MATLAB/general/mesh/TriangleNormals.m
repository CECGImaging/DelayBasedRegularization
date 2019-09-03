function [normals,unormals,areas]=TriangleNormals(nodes,elements)
% function [normals,unormals,areas]=TriangleNormals(nodes,elements)
%unormals = unit normals
%normals = normals with length 2A
%areas = areas of the triangles (A)

p1=nodes(elements(:,1),:);
p2=nodes(elements(:,2),:);
p3=nodes(elements(:,3),:);
normals=cross(p2-p1,p3-p1);
areas=sqrt(sum(normals.*normals,2))/2;
unormals=normals./(2*areas*[1 1 1]);

function R=cross(R1,R2)
R(:,1)=R1(:,2).*R2(:,3)-R1(:,3).*R2(:,2);
R(:,2)=R1(:,3).*R2(:,1)-R1(:,1).*R2(:,3);
R(:,3)=R1(:,1).*R2(:,2)-R1(:,2).*R2(:,1);
