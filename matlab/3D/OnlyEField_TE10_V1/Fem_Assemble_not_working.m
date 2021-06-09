% --------------------------------------------------------------
% Assemble the global matrices
% --------------------------------------------------------------
function [KMtx, BMtx, bMtx] = ...
    Fem_Assemble(no2xyz, el2no, el2ma, ma2er, ma2si, port_fac2no_list, k0, gamma, k_z10, a, E0, direction)

% Arguments:
%   no2xyz = coordinates of the nodes
%   el2no = nodes of the tetrahedrons
%   el2ma = material of the tetrahedrons
%   ma2er = relative permittivity of the materials
%   ma2si = conductivity of the materials
% Returns:

% need  to incorperate to be more flexible


%global ed2noLoc fa2noLoc
ed2noLoc = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4]';

% Global number of entities

%%%%%%%%%%%%%%%%%%%%%%%% Asemble K matrix ################

elNumGlo = size(el2no,2); % number of thetras
edNumGlo = ElementDatabase_Cardinal('edges'); % Number of edges / lines

% Incremental of data snippets for each element in idxRes_EE, idxRes_FE, idxRes_FF
incRes_EE = 6*6; % edge * edge
% Initializing. Current location based on step, starts in 1 instead of 0 due to matlab
% properties
idxRes_EE = 1;
irRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for row index 
icRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for column index 
mKRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for data

tic
noTmp = zeros([size(ed2noLoc,1),elNumGlo*size(ed2noLoc,2)]);
noTmpTmp = zeros(size(ed2noLoc));
pointer = 1;
for elIdx = 1:elNumGlo
  noTmpTmp(:) = el2no(ed2noLoc(:),elIdx);
  noTmp(:,pointer:pointer+5) = noTmpTmp;
  pointer = pointer+size(ed2noLoc,2);
end
eiVtr_list = ElementDatabase_Get('edges', noTmp) ; 

pointer = 1;

for elIdx = 1:elNumGlo % goes throug the amount of thetras
    
    no = el2no(:,elIdx); % current nodes (corners / points)
    xyz = no2xyz(:,no); % coordinates for the nodes (4 nodes)
    
    [KElMtx_EE] = ...
        Fem_Cmp_Vol_ElMtx(xyz, ma2er{el2ma(elIdx)}, ma2si{el2ma(elIdx)},k0);
    %for edges
    noTmp = zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
    noTmp(:) = el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
    esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
    %eiVtr = ElementDatabase_Get('edges', noTmp); % each edge id assosiated with the specific thetra (6 edges)
    
    
    % edges only
    irTmp_EE = eiVtr'*ones(size(eiVtr)); % row index for eMtx only based on edge id
    icTmp_EE = ones(size(eiVtr'))*eiVtr; % column index for eMtx only based on edge id
    isTmp_EE = esVtr'*esVtr; % The direction adjusted for the eElMtx_EE and sElMtx_EE matrix

    
    irRes_EE(idxRes_EE + (1:incRes_EE) - 1) = irTmp_EE(:); % includes all row indexes for all thetras for eMtx and sMtx
    icRes_EE(idxRes_EE + (1:incRes_EE) - 1) = icTmp_EE(:); % includes all column indexis for all thetras for eMtx and sMtx
    mKRes_EE(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*KElMtx_EE(:); % resulting data for eMtx

    
    idxRes_EE = idxRes_EE + incRes_EE; % incerment pointer for edge edge

end

KMtx = sparse(irRes_EE, icRes_EE, mKRes_EE, edNumGlo, edNumGlo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Asemble the B matrix, impedance boundery for all ports %%%%%%%%%%%

% Incremental of data snippets for each element in idxRes_EE, idxRes_FE, idxRes_FF
incRes_EE = 3*3; % edge * edge
% Initializing. Current location based on step, starts in 1 instead of 0 due to matlab
% properties
idxRes_EE = 1;

elNum = 0;
for i = 1:length(port_fac2no_list)
    elNum = elNum + size(port_fac2no_list{i},2);
end

irRes_EE = zeros(incRes_EE*elNum,1); % pre allocation for row index 
icRes_EE = zeros(incRes_EE*elNum,1); % pre allocation for column index 
mBRes_EE = zeros(incRes_EE*elNum,1); % pre allocation for data

for port_num = 1:length(port_fac2no_list)
    for no = port_fac2no_list{port_num} % goes throug the amount of thetras

        %no = el2no(:,elIdx); % current nodes (corners / points)
        xyz = no2xyz(:,no); % coordinates for the nodes (4 nodes)

        [BElMtx_EE] = ...
            Fem_Cmp_Surface_ElMtx(xyz,gamma,direction(port_num));
        %for edges
        noTmp = zeros(2,3); %zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
        noTmp(:) = no([[1,2];[2,3];[3,1]]');%el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
        esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
        eiVtr = ElementDatabase_Get('edges', noTmp); % each edge id assosiated with the specific thetra (6 edges)

        % edges only
        irTmp_EE = eiVtr'*ones(size(eiVtr)); % row index for BMtx only based on edge id
        icTmp_EE = ones(size(eiVtr'))*eiVtr; % column index for BMtx only based on edge id
        isTmp_EE = esVtr'*esVtr; % The direction adjusted for the eElMtx_EE and sElMtx_EE matrix


        irRes_EE(idxRes_EE + (1:incRes_EE) - 1) = irTmp_EE(:); % includes all row indexes for all thetras for eMtx and sMtx
        icRes_EE(idxRes_EE + (1:incRes_EE) - 1) = icTmp_EE(:); % includes all column indexis for all thetras for eMtx and sMtx
        mBRes_EE(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*BElMtx_EE(:); % resulting data for eMtx

        idxRes_EE = idxRes_EE + incRes_EE; % incerment pointer for edge edge

    end
end

BMtx = sparse(irRes_EE, icRes_EE, mBRes_EE, edNumGlo, edNumGlo);

% %%%%%%%%%%%%%%%%%%%% small b active port %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        Fem_Cmp_Surface_active_ElMtx(xyz,k_z10,a,E0,direction(1));
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
bMtx = sparse(irRes_EE, icRes_EE, mbRes_EE, edNumGlo, 1);