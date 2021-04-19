% --------------------------------------------------------------
% Assemble the global matrices
% --------------------------------------------------------------
function [KMtx, BMtx] = ...
    Fem_Assemble(no2xyz, el2no, el2ma, ma2er, ma2si, fac2no_port1, fac2no_port2)

% Arguments:
%   no2xyz = coordinates of the nodes
%   el2no = nodes of the tetrahedrons
%   el2ma = material of the tetrahedrons
%   ma2er = relative permittivity of the materials
%   ma2si = conductivity of the materials
% Returns:

% need  to incorperate to be more flexible
a = 0.2;
c0 = 299792458;
f = 2*10^9; % 2 ghz
k0 = (f/c0)^2 * 4*pi^2;

% only needed for small b
k_z10 = sqrt(k0^2-(pi/a)^2);
gamma = 1j*k_z10;



global ed2noLoc fa2noLoc

% Global number of entities
elNumGlo = size(el2no,2); % number of thetras
edNumGlo = ElementDatabase_Cardinal('edges'); % Number of edges / lines

% Incremental of data snippets for each element in idxRes_EE, idxRes_FE, idxRes_FF
incRes_EE = 6*6; % edge * edge
% Initializing. Current location based on step, starts in 1 instead of 0 due to matlab
% properties
idxRes_EE = 1;

irRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for row index to spase matrix eMtx and sMtx
icRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for column index to spase matrix eMtx and sMtx

meRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for data in eMtx

for elIdx = 1:elNumGlo % goes throug the amount of thetras
    
    no = el2no(:,elIdx); % current nodes (corners / points)
    xyz = no2xyz(:,no); % coordinates for the nodes (4 nodes)
    
    [KElMtx_EE] = ...
        Fem_Cmp_Vol_ElMtx(xyz, ma2er{el2ma(elIdx)}, ma2si{el2ma(elIdx)},k0);
    %for edges
    noTmp = zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
    noTmp(:) = el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
    esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
    eiVtr = ElementDatabase_Get('edges', noTmp); % each edge id assosiated with the specific thetra (6 edges)

    % edges only
    irTmp_EE = eiVtr'*ones(size(eiVtr)); % row index for eMtx only based on edge id
    icTmp_EE = ones(size(eiVtr'))*eiVtr; % column index for eMtx only based on edge id
    isTmp_EE = esVtr'*esVtr; % The direction adjusted for the eElMtx_EE and sElMtx_EE matrix

    
    irRes_EE(idxRes_EE + (1:incRes_EE) - 1) = irTmp_EE(:); % includes all row indexes for all thetras for eMtx and sMtx
    icRes_EE(idxRes_EE + (1:incRes_EE) - 1) = icTmp_EE(:); % includes all column indexis for all thetras for eMtx and sMtx
    meRes_EE(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*KElMtx_EE(:); % resulting data for eMtx

    
    idxRes_EE = idxRes_EE + incRes_EE; % incerment pointer for edge edge

end

KMtx = sparse(irRes_EE, icRes_EE, meRes_EE, edNumGlo, edNumGlo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to change this later so that it is adaptive depending on how manny
% ports are assigned
% edIdx_port1 = ElementDatabase_Get('edges', ed2no_port1);
% edIdx_port2 = ElementDatabase_Get('edges', ed2no_port2);

% Incremental of data snippets for each element in idxRes_EE, idxRes_FE, idxRes_FF
incRes_EE = 3*3; % edge * edge

% Initializing. Current location based on step, starts in 1 instead of 0 due to matlab
% properties
idxRes_EE = 1;

irRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for row index to spase matrix eMtx and sMtx
icRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for column index to spase matrix eMtx and sMtx


meRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for data in eMtx

for no = [fac2no_port1, fac2no_port2] % goes throug the amount of thetras
    
    %no = el2no(:,elIdx); % current nodes (corners / points)
    xyz = no2xyz(:,no); % coordinates for the nodes (4 nodes)
    
    [bElMtx_EE,BElMtx_EE] = ...
        Fem_Cmp_Surface_ElMtx(xyz,gamma);
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
    mbRes_EE(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*bElMtx_EE(:); % resulting data for eMtx

    
    idxRes_EE = idxRes_EE + incRes_EE; % incerment pointer for edge edge

end


BMtx = sparse(irRes_EE, icRes_EE, msRes_EE, edNumGlo, edNumGlo);
%bMtx =


% %%%%%%%%%%%%%%%%%%%% small b %%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% incRes_EE = 3*3; % edge * edge
% 
% % Initializing. Current location based on step, starts in 1 instead of 0 due to matlab
% % properties
% idxRes_EE = 1;
% 
% irRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for row index to spase matrix eMtx and sMtx
% icRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for column index to spase matrix eMtx and sMtx
% 
% 
% meRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for data in eMtx
% 
% for no = [fac2no_port1, fac2no_port2] % goes throug the amount of thetras
%     
%     %no = el2no(:,elIdx); % current nodes (corners / points)
%     xyz = no2xyz(:,no); % coordinates for the nodes (4 nodes)
%     
%     [KElMtx_EE] = ...
%         Fem_Cmp_Surface_ElMtx(xyz,gamma);
%     %for edges
%     noTmp = zeros(2,3); %zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
%     noTmp(:) = no([[1,2];[2,3];[3,1]]');%el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
%     esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
%     eiVtr = ElementDatabase_Get('edges', noTmp); % each edge id assosiated with the specific thetra (6 edges)
% 
%     % edges only
%     irTmp_EE = eiVtr'*ones(size(eiVtr)); % row index for BMtx only based on edge id
%     icTmp_EE = ones(size(eiVtr'))*eiVtr; % column index for BMtx only based on edge id
%     isTmp_EE = esVtr'*esVtr; % The direction adjusted for the eElMtx_EE and sElMtx_EE matrix
% 
%     
%     irRes_EE(idxRes_EE + (1:incRes_EE) - 1) = irTmp_EE(:); % includes all row indexes for all thetras for eMtx and sMtx
%     icRes_EE(idxRes_EE + (1:incRes_EE) - 1) = icTmp_EE(:); % includes all column indexis for all thetras for eMtx and sMtx
%     meRes_EE(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*KElMtx_EE(:); % resulting data for eMtx
% 
%     
%     idxRes_EE = idxRes_EE + incRes_EE; % incerment pointer for edge edge
% 
% end
% 
% 
% BMtx = sparse(irRes_EE, icRes_EE, msRes_EE, edNumGlo, edNumGlo);