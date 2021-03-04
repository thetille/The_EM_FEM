% --------------------------------------------------------------
% Compute the stiffness and mass matrix for edge elements on
% a triangular grid
% --------------------------------------------------------------
function [M, S, el2ed] = edgeFEM2D(no2xy, el2no)
% Arguments:
% no2xy = x- and y-coordinates of the nodes
% el2no = node indices of the triangles
% Returns:
% M = Mass matrix
% S = Stiffness matrix
% el2ed = a table that contain the three edge numbers related
% to each element
% Sort the nodes of each element
el2no = sort(el2no);
% Assign a number to each edge in the grid and create el2ed
n1 = el2no([1 1 2],:);
n2 = el2no([2 3 3],:);
[ed2no,trash,el2ed] = unique([n1(:) n2(:)],'rows');
el2ed = reshape(el2ed,3,size(el2no,2));
% Compute det(Jˆe), grad phi_2 and grad phi_3
e1 = no2xy(:,el2no(2,:)) - no2xy(:,el2no(1,:)); % 1st edge in
% all elements
e2 = no2xy(:,el2no(3,:)) - no2xy(:,el2no(1,:)); % 2nd edge in
% all elements
detJ = e1(1,:).*e2(2,:) - e1(2,:).*e2(1,:); % det(Jˆe) for
% all elements
g2 = [+e2(2,:)./detJ; -e2(1,:)./detJ]; % grad phi_2
g3 = [-e1(2,:)./detJ; +e1(1,:)./detJ]; % grad phi_3
% Define element shape independent matrices
m22 = [+3 +1 -1; +1 +1 -1; -1 -1 +1] / 12;
m23 = [+3 +3 +1; +3 +3 -1; +1 -1 -1] / 12;
m33 = [+1 +1 +1; +1 +3 +1; +1 +1 +1] / 12;
s00 = [+2 -2 +2; -2 +2 -2; +2 -2 +2];
% Compute local matrices and indices for all elements
mloc = m22(:) * (abs(detJ).*sum(g2.*g2)) + ...
m23(:) * (abs(detJ).*sum(g2.*g3)) + ...
m33(:) * (abs(detJ).*sum(g3.*g3));
sloc = s00(:) * abs(1./detJ);
rows = el2ed([1 2 3 1 2 3 1 2 3],:);
cols = el2ed([1 1 1 2 2 2 3 3 3],:);
% Assemble.
S = sparse(rows,cols,sloc);
M = sparse(rows,cols,mloc);