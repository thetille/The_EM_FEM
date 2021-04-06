% --------------------------------------------------------------
% Projection from the curl-conforming linear elements on 
% tetrahedrons to the space of piecewise linear node based 
% elements. The projection matrices are used as
%
%    exVtr = proj_ed2noMtx.xc * eVtr
%    eyVtr = proj_ed2noMtx.yc * eVtr
%    ezVtr = proj_ed2noMtx.zc * eVtr
%
% where 'eVtr' is the degrees of freedom associated with 
% the edge elements. Further, 'exVtr' is the degrees of freedom 
% associated with the nodes that describe the x-component of 
% the electric field. Similarily, 'eyVtr' and 'ezVtr' correspond 
% to the y- and z-components. For the magnetic flux density, we have
%
%      bxVtr = proj_fa2noMtx.xc * bVtr
%      byVtr = proj_fa2noMtx.yc * bVtr
%      bzVtr = proj_fa2noMtx.zc * bVtr
% --------------------------------------------------------------
function [proj_ed2noMtx, proj_fa2noMtx] = ...
    ProjSol2Nodes_Assemble(no2xyz, el2no)

% Arguments:
%   no2xyz = coordinates of the nodes
%   el2no = nodes of the elements
% Returns:
%   proj_ed2noMtx = projection matrix from edge elements to nodes
%   proj_fa2noMtx = projection matrix from face elements to nodes

global ed2noLoc fa2noLoc

% Global number of entities
elNumGlo = size(el2no,2);
noNumGlo = size(no2xyz,2);
edNumGlo = ElementDatabase_Cardinal('edges');
faNumGlo = ElementDatabase_Cardinal('faces');

% Incremental steps for each element
incRes_NN = 4*4;
incRes_NE = 4*6;
incRes_NF = 4*4;

% Initializing global counters
idxRes_NN = 1;
idxRes_NE = 1;
idxRes_NF = 1;

irRes_NN = zeros(incRes_NN*elNumGlo,1);
icRes_NN = zeros(incRes_NN*elNumGlo,1);
msRes_NN = zeros(incRes_NN*elNumGlo,1);

irRes_NE = zeros(incRes_NE*elNumGlo,1);
icRes_NE = zeros(incRes_NE*elNumGlo,1);
mxRes_NE = zeros(incRes_NE*elNumGlo,1);
myRes_NE = zeros(incRes_NE*elNumGlo,1);
mzRes_NE = zeros(incRes_NE*elNumGlo,1);

irRes_NF = zeros(incRes_NF*elNumGlo,1);
icRes_NF = zeros(incRes_NF*elNumGlo,1);
mxRes_NF = zeros(incRes_NF*elNumGlo,1);
myRes_NF = zeros(incRes_NF*elNumGlo,1);
mzRes_NF = zeros(incRes_NF*elNumGlo,1);

% Computing the contributions to the mass and stiffness matrices.
for elIdx = 1:elNumGlo
    
    no = el2no(:,elIdx);
    xyz = no2xyz(:,no);
    
    [iElMtx_NN, iElMtx_NE, iElMtx_NF] = ProjSol2Nodes_CmpElMtx(xyz);
    
    nsVtr = ones(size(no'));
    niVtr = no';    
    
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
    
    irTmp_NN = niVtr'*ones(size(niVtr));
    icTmp_NN = ones(size(niVtr'))*niVtr;
    
    irTmp_NE = niVtr'*ones(size(eiVtr));
    icTmp_NE = ones(size(niVtr'))*eiVtr;
    isTmp_NE = nsVtr'*esVtr;
    
    irTmp_NF = niVtr'*ones(size(fiVtr));
    icTmp_NF = ones(size(niVtr'))*fiVtr;
    isTmp_NF = nsVtr'*fsVtr;
    
    irRes_NN(idxRes_NN + (1:incRes_NN) - 1) = irTmp_NN(:);
    icRes_NN(idxRes_NN + (1:incRes_NN) - 1) = icTmp_NN(:);
    msRes_NN(idxRes_NN + (1:incRes_NN) - 1) = iElMtx_NN(:);
    
    irRes_NE(idxRes_NE + (1:incRes_NE) - 1) = irTmp_NE(:);
    icRes_NE(idxRes_NE + (1:incRes_NE) - 1) = icTmp_NE(:);
    mxRes_NE(idxRes_NE + (1:incRes_NE) - 1) = isTmp_NE(:).*iElMtx_NE.x(:);
    myRes_NE(idxRes_NE + (1:incRes_NE) - 1) = isTmp_NE(:).*iElMtx_NE.y(:);
    mzRes_NE(idxRes_NE + (1:incRes_NE) - 1) = isTmp_NE(:).*iElMtx_NE.z(:);

    irRes_NF(idxRes_NF + (1:incRes_NF) - 1) = irTmp_NF(:);
    icRes_NF(idxRes_NF + (1:incRes_NF) - 1) = icTmp_NF(:);
    mxRes_NF(idxRes_NF + (1:incRes_NF) - 1) = isTmp_NF(:).*iElMtx_NF.x(:);
    myRes_NF(idxRes_NF + (1:incRes_NF) - 1) = isTmp_NF(:).*iElMtx_NF.y(:);
    mzRes_NF(idxRes_NF + (1:incRes_NF) - 1) = isTmp_NF(:).*iElMtx_NF.z(:);

    idxRes_NN = idxRes_NN + incRes_NN;
    idxRes_NE = idxRes_NE + incRes_NE;
    idxRes_NF = idxRes_NF + incRes_NF;
    
end

msMtx_NN = sparse(irRes_NN, icRes_NN, msRes_NN, noNumGlo, noNumGlo);

mxMtx_NE = sparse(irRes_NE, icRes_NE, mxRes_NE, noNumGlo, edNumGlo);
myMtx_NE = sparse(irRes_NE, icRes_NE, myRes_NE, noNumGlo, edNumGlo);
mzMtx_NE = sparse(irRes_NE, icRes_NE, mzRes_NE, noNumGlo, edNumGlo);

mxMtx_NF = sparse(irRes_NF, icRes_NF, mxRes_NF, noNumGlo, faNumGlo);
myMtx_NF = sparse(irRes_NF, icRes_NF, myRes_NF, noNumGlo, faNumGlo);
mzMtx_NF = sparse(irRes_NF, icRes_NF, mzRes_NF, noNumGlo, faNumGlo);

isMtx = diag(1./diag(msMtx_NN));

proj_ed2noMtx.xc = isMtx*mxMtx_NE;
proj_ed2noMtx.yc = isMtx*myMtx_NE;
proj_ed2noMtx.zc = isMtx*mzMtx_NE;

proj_fa2noMtx.xc = isMtx*mxMtx_NF;
proj_fa2noMtx.yc = isMtx*myMtx_NF;
proj_fa2noMtx.zc = isMtx*mzMtx_NF;


