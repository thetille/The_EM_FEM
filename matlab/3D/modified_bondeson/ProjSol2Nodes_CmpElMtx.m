% --------------------------------------------------------------
% Compute projection matrices that project curl-conforming
% and divergence-conforming field solutions to the nodes
% of the tetrahedrons
% --------------------------------------------------------------
function [iElMtx_NN, iElMtx_NE, iElMtx_NF] = ...
    ProjSol2Nodes_CmpElMtx(xyz)

% Argument:
%   xyz = the coordinates of the nodes of the element
%   xyz(1,:) = x-coordinate of all nodes
%   xyz(2,:) = y-coordinate of all nodes
%   xyz(3,:) = z-coordinate of all nodes
%   xyz(:,n) = xyz-coordinates of the n:th node
% Returns:
%   iElMtx_NN = element projection matrix (phi_i phi_j)
%   iElMtx_NE = element projection matrix (phi_i N_j)
%   iElMtx_NF = element projection matrix (phi_i M_j)

% Quadrature rule
q2r = [[5.854101966249685e-01, ...
    1.381966011250105e-01, ...
    1.381966011250105e-01]; ...
    [1.381966011250105e-01, ...
    5.854101966249685e-01, ...
    1.381966011250105e-01]; ...
    [1.381966011250105e-01, ...
    1.381966011250105e-01, ...
    5.854101966249685e-01]; ...
    [1.381966011250105e-01, ...
    1.381966011250105e-01, ...
    1.381966011250105e-01]]';
q2w = [4.166666666666666e-02; ...
    4.166666666666666e-02; ...
    4.166666666666666e-02; ...
    4.166666666666666e-02]';

q2r_lump = [[0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00]; ...
    [1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00]; ...
    [0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00]; ...
    [0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00]]';
q2w_lump = [4.166666666666666e-02; ...
    4.166666666666666e-02; ...
    4.166666666666666e-02; ...
    4.166666666666666e-02]';

% H(grad) basis functions
up{1} = 1 - q2r(1,:) - q2r(2,:) - q2r(3,:);
up{2} = q2r(1,:);
up{3} = q2r(2,:);
up{4} = q2r(3,:);

up_lump{1} = 1 - q2r_lump(1,:) - q2r_lump(2,:) - q2r_lump(3,:);
up_lump{2} = q2r_lump(1,:);
up_lump{3} = q2r_lump(2,:);
up_lump{4} = q2r_lump(3,:);

% Gradient of H(grad) basis functions
ug{1} = [-1 -1 -1]';
ug{2} = [+1 0 0]';
ug{3} = [0 +1 0]';
ug{4} = [0 0 +1]';

% H(div) basis functions
uim{1} = 2*(cross(ug{2},ug{1})*up{3} + cross(ug{1},ug{3})*up{2} + cross(ug{3},ug{2})*up{1});
uim{2} = 2*(cross(ug{2},ug{4})*up{1} + cross(ug{4},ug{1})*up{2} + cross(ug{1},ug{2})*up{4});
uim{3} = 2*(cross(ug{3},ug{4})*up{2} + cross(ug{4},ug{2})*up{3} + cross(ug{2},ug{3})*up{4});
uim{4} = 2*(cross(ug{1},ug{4})*up{3} + cross(ug{4},ug{3})*up{1} + cross(ug{3},ug{1})*up{4});

% H(curl) basis function
uin{1} = ug{2}*up{1} - ug{1}*up{2};
uin{2} = ug{3}*up{2} - ug{2}*up{3};
uin{3} = ug{1}*up{3} - ug{3}*up{1};
uin{4} = ug{4}*up{1} - ug{1}*up{4};
uin{5} = ug{4}*up{2} - ug{2}*up{4};
uin{6} = ug{4}*up{3} - ug{3}*up{4};

% Curl of H(curl) basis functions
ouTmp = ones(size(q2w));
ucn{1} = 2*cross(ug{1},ug{2})*ouTmp;
ucn{2} = 2*cross(ug{2},ug{3})*ouTmp;
ucn{3} = 2*cross(ug{3},ug{1})*ouTmp;
ucn{4} = 2*cross(ug{1},ug{4})*ouTmp;
ucn{5} = 2*cross(ug{2},ug{4})*ouTmp;
ucn{6} = 2*cross(ug{3},ug{4})*ouTmp;

% Jacobian
jac = zeros(3);
for iIdx = 1:4
    jac = jac ...
        + [xyz(1,iIdx)*ug{iIdx}, ...
        xyz(2,iIdx)*ug{iIdx}, ...
        xyz(3,iIdx)*ug{iIdx}];
end

% Mappings
det_jac = det(jac);
map_ccs = inv(jac);      % mapping for curl-conforming space
map_dcs = jac'/det_jac;  % mapping for div-conforming space
for iIdx = 1:6
    gin{iIdx} = map_ccs*uin{iIdx};
    gcn{iIdx} = map_dcs*ucn{iIdx};
end
for iIdx = 1:4
    gim{iIdx} = map_dcs*uim{iIdx};
end

% Evaluation of element matrix: phi_i phi_j
iElMtx_NN = zeros(4);
for iIdx = 1:4
    for jIdx = 1:4
        bsiTmp = up_lump{iIdx};
        bsjTmp = up_lump{jIdx};
        iElMtx_NN(iIdx,jIdx) = (bsiTmp.*bsjTmp) * q2w' * det_jac;
    end
end

% Evaluation of element matrix: phi_i (x * Nj)
iElMtx_NE.x = zeros(4,6);
iElMtx_NE.y = zeros(4,6);
iElMtx_NE.z = zeros(4,6);
for iIdx = 1:4
    for jIdx = 1:6
        bsiTmp = up{iIdx};
        bxjTmp = gin{jIdx}(1,:);
        byjTmp = gin{jIdx}(2,:);
        bzjTmp = gin{jIdx}(3,:);
        iElMtx_NE.x(iIdx,jIdx) = (bsiTmp.*bxjTmp) * q2w' * det_jac;
        iElMtx_NE.y(iIdx,jIdx) = (bsiTmp.*byjTmp) * q2w' * det_jac;
        iElMtx_NE.z(iIdx,jIdx) = (bsiTmp.*bzjTmp) * q2w' * det_jac;
    end
end

% Evaluation of element matrix: phi_i (x * Mj)
iElMtx_NF.x = zeros(4,4);
iElMtx_NF.y = zeros(4,4);
iElMtx_NF.z = zeros(4,4);
for iIdx = 1:4
    for jIdx = 1:4
        bsiTmp = up{iIdx};
        bxjTmp = gim{jIdx}(1,:);
        byjTmp = gim{jIdx}(2,:);
        bzjTmp = gim{jIdx}(3,:);
        iElMtx_NF.x(iIdx,jIdx) = (bsiTmp.*bxjTmp) * q2w' * det_jac;
        iElMtx_NF.y(iIdx,jIdx) = (bsiTmp.*byjTmp) * q2w' * det_jac;
        iElMtx_NF.z(iIdx,jIdx) = (bsiTmp.*bzjTmp) * q2w' * det_jac;
    end
end