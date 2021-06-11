% --------------------------------------------------------------
% Assemble the global matrices
% --------------------------------------------------------------
function [KMtx1,KMtx2] = ...
    Assemble_K(no2xyz, el2no, el2ma, ma2er, ma2mu)


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
mKRes_EE1 = zeros(incRes_EE*elNumGlo,1); % pre allocation for data
mKRes_EE2 = zeros(incRes_EE*elNumGlo,1); % pre allocation for data


noTmp = zeros([size(ed2noLoc,1),elNumGlo*size(ed2noLoc,2)]);
noTmpTmp = zeros(size(ed2noLoc));
pointer = 1;
for elIdx = 1:elNumGlo
  noTmpTmp(:) = el2no(ed2noLoc(:),elIdx);
  noTmp(:,pointer:pointer+5) = noTmpTmp;
  pointer = pointer+size(ed2noLoc,2);
end
eiVtr_list = ElementDatabase_Get2('edges', noTmp) ; 

pointer = 1;

for elIdx = 1:elNumGlo % goes throug the amount of thetras
    
    no = el2no(:,elIdx); % current nodes (corners / points)
    xyz = no2xyz(:,no); % coordinates for the nodes (4 nodes)
    
    [KElMtx_EE1,KElMtx_EE2] = ...
        Fem_Cmp_Vol_ElMtx(xyz, ma2er{el2ma(elIdx)}, ma2mu{el2ma(elIdx)});
    %for edges
    noTmp = zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
    noTmp(:) = el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
    esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
    %eiVtr = ElementDatabase_Get('edges', noTmp); % each edge id assosiated with the specific thetra (6 edges)
    eiVtr = eiVtr_list(pointer:pointer+5);
    
    
    % edges only
    irTmp_EE = eiVtr'*ones(size(eiVtr)); % row index for eMtx only based on edge id
    icTmp_EE = ones(size(eiVtr'))*eiVtr; % column index for eMtx only based on edge id
    isTmp_EE = esVtr'*esVtr; % The direction adjusted for the eElMtx_EE and sElMtx_EE matrix

    
    irRes_EE(idxRes_EE + (1:incRes_EE) - 1) = irTmp_EE(:); % includes all row indexes for all thetras for eMtx and sMtx
    icRes_EE(idxRes_EE + (1:incRes_EE) - 1) = icTmp_EE(:); % includes all column indexis for all thetras for eMtx and sMtx
    mKRes_EE1(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*KElMtx_EE1(:); % resulting data for eMtx
    mKRes_EE2(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*KElMtx_EE2(:);

    
    idxRes_EE = idxRes_EE + incRes_EE; % incerment pointer for edge edge
    pointer = pointer+6;
end

KMtx1 = sparse(irRes_EE, icRes_EE, mKRes_EE1, edNumGlo, edNumGlo);
KMtx2 = sparse(irRes_EE, icRes_EE, mKRes_EE2, edNumGlo, edNumGlo);