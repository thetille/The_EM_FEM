% --------------------------------------------------------------
% Initialize the FE-solver by numbering of local and global
% edges and faces.
% --------------------------------------------------------------
function Fem_Init(no2xyz, ed2no_all, fa2no_all)

% Arguments:
%    no2xyz = coordinates of the nodes
%    ed2no_all = nodes of all edges, lines
%    fa2no_all = nodes of all faces, triangles
% Returns:
%    -

global ed2noLoc fa2noLoc % global vairbles in the remaining code i.e no need to have a return

% Setting up the edge information.
ed2noLoc = ...
    [1 2; 2 3; 3 1; 1 4; 2 4; 3 4]'; % base lines
fa2noLoc = ...
    [3 2 1; 1 2 4; 2 3 4; 3 1 4]';  % base triangles

% Number the edges
ElementDatabase_Init('edges', size(no2xyz,2)*[1 1])
ElementDatabase_Set('edges', ed2no_all)

% Number the faces
ElementDatabase_Init('faces', size(no2xyz,2)*[1 1 1])
ElementDatabase_Set('faces', fa2no_all)