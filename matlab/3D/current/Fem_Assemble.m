% --------------------------------------------------------------
% Assemble the global matrices
% --------------------------------------------------------------
function [eMtx, sMtx, cMtx, uMtx] = ...
    Fem_Assemble(no2xyz, el2no, el2ma, ma2er, ma2si)

% Arguments:
%   no2xyz = coordinates of the nodes
%   el2no = nodes of the tetrahedrons
%   el2ma = material of the tetrahedrons
%   ma2er = relative permittivity of the materials
%   ma2si = conductivity of the materials
% Returns:
%   eMtx = mass matrix with permittivity coefficient
%   sMtx = mass matrix with conductivity coefficient
%   cMtx = curl matrix
%   uMtx = mass matrix with unity coefficient

global ed2noLoc fa2noLoc

% Global number of entities
elNumGlo = size(el2no,2); % number of thetras
edNumGlo = ElementDatabase_Cardinal('edges'); % Number of edges / lines
faNumGlo = ElementDatabase_Cardinal('faces'); % number of faces / triangles

% Incremental of data snippets for each element in idxRes_EE, idxRes_FE, idxRes_FF
incRes_EE = 6*6; % edge * edge
incRes_FE = 4*6; % face * edge 
incRes_FF = 4*4; % face * face

% Initializing. Current location based on step, starts in 1 instead of 0 due to matlab
% properties
idxRes_EE = 1;
idxRes_FE = 1;
idxRes_FF = 1;

irRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for row index to spase matrix eMtx and sMtx
icRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for column index to spase matrix eMtx and sMtx

irRes_FE = zeros(incRes_FE*elNumGlo,1); % pre allocation for row index to spase matrix cMtx
icRes_FE = zeros(incRes_FE*elNumGlo,1); % pre allocation for column index to spase matrix cMtx

irRes_FF = zeros(incRes_FF*elNumGlo,1); % pre allocation for row index to spase matrix uMtx
icRes_FF = zeros(incRes_FF*elNumGlo,1); % pre allocation for column index to spase matrix uMtx

meRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for data in eMtx
msRes_EE = zeros(incRes_EE*elNumGlo,1); % pre allocation for data in sMtx
mcRes_FE = zeros(incRes_FE*elNumGlo,1); % pre allocation for data in cMtx
muRes_FF = zeros(incRes_FF*elNumGlo,1); % pre allocation for data in uMtx

% Computing the contributions to the mass and
% stiffness matrices.
for elIdx = 1:elNumGlo % goes throug the amount of thetras
    
    no = el2no(:,elIdx); % current nodes (corners / points)
    xyz = no2xyz(:,no); % coordinates for the nodes (4 nodes)
    
    [eElMtx_EE, sElMtx_EE, cElMtx_FE, uElMtx_FF] = ...
        Fem_CmpElMtx(xyz, ma2er{el2ma(elIdx)}, ma2si{el2ma(elIdx)});
    %for edges
    noTmp = zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
    noTmp(:) = el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
    esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
    eiVtr = ElementDatabase_Get('edges', noTmp); % each edge assosiated with the specific thetra
    
    %for faces
    noTmp = zeros(size(fa2noLoc)); % temp var of size of the amount of initial faces
    noTmp(:) = el2no(fa2noLoc(:),elIdx); % temp var with initial nodes of each face
    fsVtr = 2*(... % each edge gets a normal vector, initial direction is based on index only
        ((noTmp(1,:) < noTmp(2,:)) & (noTmp(2,:) < noTmp(3,:))) | ...
        ((noTmp(2,:) < noTmp(3,:)) & (noTmp(3,:) < noTmp(1,:))) | ...
        ((noTmp(3,:) < noTmp(1,:)) & (noTmp(1,:) < noTmp(2,:))) ...
        ) - 1;
    fiVtr = ElementDatabase_Get('faces', noTmp); % gets ids for faces assosiated with the specific thetra
        
    irTmp_EE = eiVtr'*ones(size(eiVtr));
    icTmp_EE = ones(size(eiVtr'))*eiVtr;
    isTmp_EE = esVtr'*esVtr;
    
    irTmp_FE = fiVtr'*ones(size(eiVtr));
    icTmp_FE = ones(size(fiVtr'))*eiVtr;
    isTmp_FE = fsVtr'*esVtr;
    
    irTmp_FF = fiVtr'*ones(size(fiVtr));
    icTmp_FF = ones(size(fiVtr'))*fiVtr;
    isTmp_FF = fsVtr'*fsVtr;
    
    irRes_EE(idxRes_EE + (1:incRes_EE) - 1) = irTmp_EE(:);
    icRes_EE(idxRes_EE + (1:incRes_EE) - 1) = icTmp_EE(:);
    meRes_EE(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*eElMtx_EE(:);
    msRes_EE(idxRes_EE + (1:incRes_EE) - 1) = isTmp_EE(:).*sElMtx_EE(:);
    
    irRes_FE(idxRes_FE + (1:incRes_FE) - 1) = irTmp_FE(:);
    icRes_FE(idxRes_FE + (1:incRes_FE) - 1) = icTmp_FE(:);
    mcRes_FE(idxRes_FE + (1:incRes_FE) - 1) = isTmp_FE(:).*cElMtx_FE(:);
    
    irRes_FF(idxRes_FF + (1:incRes_FF) - 1) = irTmp_FF(:);
    icRes_FF(idxRes_FF + (1:incRes_FF) - 1) = icTmp_FF(:);
    muRes_FF(idxRes_FF + (1:incRes_FF) - 1) = isTmp_FF(:).*uElMtx_FF(:);
    
    idxRes_EE = idxRes_EE + incRes_EE;
    idxRes_FE = idxRes_FE + incRes_FE;
    idxRes_FF = idxRes_FF + incRes_FF;
    
end

eMtx = sparse(irRes_EE, icRes_EE, meRes_EE, edNumGlo, edNumGlo);
sMtx = sparse(irRes_EE, icRes_EE, msRes_EE, edNumGlo, edNumGlo);
cMtx = sparse(irRes_FE, icRes_FE, mcRes_FE, faNumGlo, edNumGlo);
uMtx = sparse(irRes_FF, icRes_FF, muRes_FF, faNumGlo, faNumGlo);