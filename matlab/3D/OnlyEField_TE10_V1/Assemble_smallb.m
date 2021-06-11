% --------------------------------------------------------------
% Assemble the global matrices
% --------------------------------------------------------------
function [bMtx] = ...
    Assemble_smallb(no2xyz, port_fac2no_list, a, direction)

%global ed2noLoc fa2noLoc
ed2noLoc = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4]';

% %%%%%%%%%%%%%%%%%%%% small b active port %%%%%%%%%%%%%%%%%%%%%%%%%%%


% Incremental of data snippets for each element in idxRes_EE, idxRes_FE, idxRes_FF
incRes_EE = 3; % edge
% Initializing. Current location based on step, starts in 1 instead of 0 due to matlab
% properties
idxRes_EE = 1;

elNum = size(port_fac2no_list{1},2);
irRes_EE = zeros(incRes_EE*elNum,1); % pre allocation for row index 
icRes_EE = ones(incRes_EE*elNum,1); % pre allocation for column index 
mbRes_EE = zeros(incRes_EE*elNum,1); % pre allocation for data

global temp1;
global temp2;

temp1 = zeros(3,3);
temp2 = zeros(3,3);

for no = port_fac2no_list{1} % goes throug the amount of thetras
    
    %no = el2no(:,elIdx); % current nodes (corners / points)
    xyz = no2xyz(:,no); % coordinates for the nodes (4 nodes)
    
    [bElMtx_EE] = ...
        Fem_Cmp_Surface_active_ElMtx(xyz,a,direction(1));
    %for edges
    noTmp = zeros(2,3); %zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
    noTmp(:) = no([[1,2];[2,3];[3,1]]');%el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
    esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
    eiVtr = ElementDatabase_Get('edges', noTmp); % each edge id assosiated with the specific thetra (6 edges)

    % edges only
    irTmp_EE = eiVtr';%*ones(size(eiVtr)); % row index for BMtx only based on edge id
    isTmp_EE = esVtr';%*esVtr; % The direction adjusted for the eElMtx_EE and sElMtx_EE matrix

    
    irRes_EE(idxRes_EE + (1:incRes_EE) - 1) = irTmp_EE(:); % includes all row indexes for all thetras for eMtx and sMtx
    mbRes_EE(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*bElMtx_EE(:); % resulting data for eMtx
    
    idxRes_EE = idxRes_EE + incRes_EE; % incerment pointer for edge edge

end

edNumGlo = ElementDatabase_Cardinal('edges'); % Number of edges / lines

bMtx = sparse(irRes_EE, icRes_EE, mbRes_EE, edNumGlo, 1);