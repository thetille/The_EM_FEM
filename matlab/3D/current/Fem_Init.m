% --------------------------------------------------------------
% Initialize the FE-solver by numbering of local and global
% edges and faces.
% --------------------------------------------------------------
function Fem_Init(no2xyz, ed2no_all, fa2no_all)

% Arguments:
%    no2xyz = coordinates of the nodes
%    ed2no_all = nodes of all edges
%    fa2no_all = nodes of all faces
% Returns:
%    -

global ed2noLoc fa2noLoc

% Setting up the edge information.
ed2noLoc = ...
    [1 2; 2 3; 3 1; 1 4; 2 4; 3 4]';
fa2noLoc = ...
    [3 2 1; 1 2 4; 2 3 4; 3 1 4]';

% Number the edges
ElementDatabase_Init('edges', size(no2xyz,2)*[1 1])
ElementDatabase_Set('edges', ed2no_all)

% Number the faces
ElementDatabase_Init('faces', size(no2xyz,2)*[1 1 1])
ElementDatabase_Set('faces', fa2no_all)