% --------------------------------------------------------------
%  Mesh stl file to tethrahedrons
% --------------------------------------------------------------
% Return values:
%    no2xyz = coordinates of the nodes
%    ed2no_all = nodes of all edges, lines
%    fa2no_all = nodes of all faces, triangles
%    el2no = nodes of the tetrahedrons
%    el2ma = material of the tetrahedrons
%    ed2no_pec = edge to node for PEC
function tethrahedron_mesh(file_name)
% import stl file
structure = importGeometry(file_name);

% plot stl file
figure
pdegplot(structure,'FaceLabels','on');

% Mesh structure




% return varibles
no2xyz
ed2no_all
fa2no_all
el2no
el2ma
ed2no_pec