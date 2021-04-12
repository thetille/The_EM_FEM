% --------------------------------------------------------------
% Compute element matrices for the tetrahedron by means of
% numerical integration on the reference element
% --------------------------------------------------------------
function [KElMtx_EE] = ...
    Fem_CmpElMtx(xyz, ma2er, ma2si, k0)

% Argument:
%   xyz = the coordinates of the nodes of the element
%   ma2er = material to permittivity
%   ma2si = material to conductivity
% Returns:   
% KElMtx_EE = 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quadrature rule
q2u = [[5.854101966249685e-01, ...
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

% H(grad) basis functions
% Node basis functions for the reference element
up{1} = 1 - q2u(1,:) - q2u(2,:) - q2u(3,:);
up{2} = q2u(1,:);
up{3} = q2u(2,:);
up{4} = q2u(3,:);

% Gradient of H(grad) basis functions, used in H(div), H(curl)
ug{1} = [-1 -1 -1]';
ug{2} = [+1 0 0]';
ug{3} = [0 +1 0]';
ug{4} = [0 0 +1]';

% % H(div) basis functions
% uim{1} = 2*(  cross(ug{2},ug{1})*up{3} ...
%             + cross(ug{1},ug{3})*up{2} ...
%             + cross(ug{3},ug{2})*up{1});
% uim{2} = 2*(  cross(ug{2},ug{4})*up{1} ...
%             + cross(ug{4},ug{1})*up{2} ...
%             + cross(ug{1},ug{2})*up{4});
% uim{3} = 2*(  cross(ug{3},ug{4})*up{2} ...
%             + cross(ug{4},ug{2})*up{3} ...
%             + cross(ug{2},ug{3})*up{4});
% uim{4} = 2*(  cross(ug{1},ug{4})*up{3} ...
%             + cross(ug{4},ug{3})*up{1} ...
%             + cross(ug{3},ug{1})*up{4});
    
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

% Physical coordinates
% Maps from refference element to physical element
q2x = zeros(3,length(q2w)); % allocate memorry
for iIdx = 1:4 % 4 edge points
    q2x = q2x + xyz(:,iIdx)*up{iIdx};
end

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



% Evaluation of element matrix: curl_Ni curl_Nj
% K_{ij}
for iIdx = 1:6
    for jIdx = 1:6
        maTmp = ones(size(q2w));
        ipTmp = maTmp.*sum(gcn{iIdx}.*gcn{jIdx});
        er = ma2er(q2x(1,:),q2x(2,:),q2x(3,:));
        endTemp = k0^2*er. * sum(gin{iIdx}. * gin{jIdx});
        KElMtx_EE(iIdx,jIdx) = (ipTmp - endTemp) * q2w' * det_jac;
    end
end