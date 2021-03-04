% --------------------------------------------------------------
% Assemble the global matrices
% --------------------------------------------------------------
function [eMtx, sMtx, cMtx, uMtx] = ...
    Fem_Assemble(no2xyz, el2no, el2ma, ma2er, ma2si)
% Arguments:
% no2xyz = coordinates of the nodes
% el2no = nodes of the tetrahedrons
% el2ma = material of the tetrahedrons
% ma2er = relative permittivity of the materials
% ma2si = conductivity of the materials
% Returns:
% eMtx = mass matrix with permittivity coefficient
% sMtx = mass matrix with conductivity coefficient
% cMtx = curl matrix
% uMtx = mass matrix with unity coefficient
global ed2noLoc fa2noLoc
% Global number of entities
elNumGlo = size(el2no,2);
edNumGlo = ElementDatabase_Cardinal('edges');
faNumGlo = ElementDatabase_Cardinal('faces');
% Incremental steps for each element
incRes_EE = 6*6;
incRes_FE = 4*6;
incRes_FF = 4*4;
% Initializing.
idxRes_EE = 1;
idxRes_FE = 1;
idxRes_FF = 1;
irRes_EE = zeros(incRes_EE*elNumGlo,1);
icRes_EE = zeros(incRes_EE*elNumGlo,1);
irRes_FE = zeros(incRes_FE*elNumGlo,1);
icRes_FE = zeros(incRes_FE*elNumGlo,1);
irRes_FF = zeros(incRes_FF*elNumGlo,1);
icRes_FF = zeros(incRes_FF*elNumGlo,1);
meRes_EE = zeros(incRes_EE*elNumGlo,1);
msRes_EE = zeros(incRes_EE*elNumGlo,1);
mcRes_FE = zeros(incRes_FE*elNumGlo,1);
mmRes_FF = zeros(incRes_FF*elNumGlo,1);
% Computing the contributions to the mass and
% stiffness matrices.
for elIdx = 1:elNumGlo
    no = el2no(:,elIdx);
    xyz = no2xyz(:,no);
    [eElMtx_EE, sElMtx_EE, cElMtx_FE, uElMtx_FF] = ...
        Fem_CmpElMtx(xyz, ma2er{el2ma(elIdx)}, ma2si{el2ma(elIdx)});
    noTmp = zeros(size(ed2noLoc));
    noTmp(:) = el2no(ed2noLoc(:),elIdx);
    esVtr = sign(noTmp(2,:)-noTmp(1,:));
    eiVtr = ElementDatabase_Get('edges', noTmp);
    noTmp = zeros(size(fa2noLoc));
    noTmp(:) = el2no(fa2noLoc(:),elIdx);
    fsVtr = 2*(...
        ((noTmp(1,:) < noTmp(2,:)) & (noTmp(2,:) < noTmp(3,:))) | ...
        ((noTmp(2,:) < noTmp(3,:)) & (noTmp(3,:) < noTmp(1,:))) | ...
        ((noTmp(3,:) < noTmp(1,:)) & (noTmp(1,:) < noTmp(2,:))) ...
        ) - 1;
    fiVtr = ElementDatabase_Get('faces', noTmp);
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