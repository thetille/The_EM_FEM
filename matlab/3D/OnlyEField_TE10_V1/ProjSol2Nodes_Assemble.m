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
function [proj_ed2noMtx] = ...
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


% Incremental steps for each element
incRes_NN = 4*4; % for nodal elements
incRes_NE = 4*6; % for node to edge elements

% Initializing global counters
idxRes_NN = 1;
idxRes_NE = 1;

irRes_NN = zeros(incRes_NN*elNumGlo,1); %placholder for row index
icRes_NN = zeros(incRes_NN*elNumGlo,1); %placeholder for colum index
msRes_NN = zeros(incRes_NN*elNumGlo,1); %placholder for results

irRes_NE = zeros(incRes_NE*elNumGlo,1); % row 
icRes_NE = zeros(incRes_NE*elNumGlo,1); % colum
mxRes_NE = zeros(incRes_NE*elNumGlo,1); % res for x
myRes_NE = zeros(incRes_NE*elNumGlo,1); % res for y
mzRes_NE = zeros(incRes_NE*elNumGlo,1); % res for z

for elIdx = 1:elNumGlo
    
    no = el2no(:,elIdx);
    xyz = no2xyz(:,no);
    
    [iElMtx_NN, iElMtx_NE] = ProjSol2Nodes_CmpElMtx(xyz);
    
    nsVtr = ones(size(no')); % direction of node elements (always one)
    niVtr = no';  %node Id
    
    noTmp = zeros(size(ed2noLoc));
    noTmp(:) = el2no(ed2noLoc(:),elIdx);
    esVtr = sign(noTmp(2,:)-noTmp(1,:)); %direction of edges
    eiVtr = ElementDatabase_Get('edges', noTmp); %Id of edges
    
    irTmp_NN = niVtr'*ones(size(niVtr)); %index row for Node to Node
    icTmp_NN = ones(size(niVtr'))*niVtr; %index colum for Node to Node
    
    irTmp_NE = niVtr'*ones(size(eiVtr)); %index row for Edges to Node
    icTmp_NE = ones(size(niVtr'))*eiVtr; %index colum for Edges to Node
    isTmp_NE = nsVtr'*esVtr; %direction for edges
    
    irRes_NN(idxRes_NN + (1:incRes_NN) - 1) = irTmp_NN(:); %save rows
    icRes_NN(idxRes_NN + (1:incRes_NN) - 1) = icTmp_NN(:); %save coulms
    msRes_NN(idxRes_NN + (1:incRes_NN) - 1) = iElMtx_NN(:); %save result for Node to Node
    
    irRes_NE(idxRes_NE + (1:incRes_NE) - 1) = irTmp_NE(:); %save rows
    icRes_NE(idxRes_NE + (1:incRes_NE) - 1) = icTmp_NE(:); %save coulms
    mxRes_NE(idxRes_NE + (1:incRes_NE) - 1) = isTmp_NE(:).*iElMtx_NE.x(:); %save x result for Node to edge
    myRes_NE(idxRes_NE + (1:incRes_NE) - 1) = isTmp_NE(:).*iElMtx_NE.y(:); %save y result for Node to edge
    mzRes_NE(idxRes_NE + (1:incRes_NE) - 1) = isTmp_NE(:).*iElMtx_NE.z(:); %save z result for Node to edge
    
    idxRes_NN = idxRes_NN + incRes_NN; %increment counter (for saving buffer node to node)
    idxRes_NE = idxRes_NE + incRes_NE; %increment counter (for saving buffer edge to node)
end

msMtx_NN = sparse(irRes_NN, icRes_NN, msRes_NN, noNumGlo, noNumGlo); %put result into sparce matrix

mxMtx_NE = sparse(irRes_NE, icRes_NE, mxRes_NE, noNumGlo, edNumGlo); %put result into sparce matrix
myMtx_NE = sparse(irRes_NE, icRes_NE, myRes_NE, noNumGlo, edNumGlo); %put result into sparce matrix
mzMtx_NE = sparse(irRes_NE, icRes_NE, mzRes_NE, noNumGlo, edNumGlo); %put result into sparce matrix

isMtx = diag(1./diag(msMtx_NN)); %invers serult of node to node

proj_ed2noMtx.xc = isMtx*mxMtx_NE; 
proj_ed2noMtx.yc = isMtx*myMtx_NE;
proj_ed2noMtx.zc = isMtx*mzMtx_NE;


