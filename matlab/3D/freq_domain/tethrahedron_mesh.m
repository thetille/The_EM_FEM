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
function [no2xyz, ed2no_all, fa2no_all, el2no, el2ma, ed2no_pec] = tethrahedron_mesh(filename)
% import mesh
[total_nodes, nodes, total_elements, triangles, tetrahedron] = import_msh_v2('waveguide_model_format2.msh')

% chanage format 
no2xyz = nodes;
ed2no_all = 






% run(filename); % Import file
% 
% no2xyz = POS';
% el2no = TETS(:,1:4)';
% ed2no_all = LINES(:,1:2)';
% fa2no_all = TRIANGKES(:,1:3)';
% 
% 
% el2ma = 0;
% ed2no_pec = 0;



% % import stl file
% model = createpde;
% structure = importGeometry(model,file_name);
% 
% % plot stl file
% figure
% pdegplot(structure,'FaceLabels','on');
% 
% % Mesh structure
% generateMesh(model);
% figure
% pdeplot3D(model);

