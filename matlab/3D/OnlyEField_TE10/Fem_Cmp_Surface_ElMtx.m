% --------------------------------------------------------------
% Compute element matrices for the port surfaces (trangles) by means of
% numerical integration on the reference element
% --------------------------------------------------------------
function [BElMtx_EE] = ...
    Fem_Cmp_Surface_ElMtx(xyz, gamma) % need to change from hardcode to adaptive
% Argument:
%   xyz = the coordinates of the nodes of the element
%   ma2er = material to permittivity
%   ma2si = material to conductivity
% Return:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% need  to incorperate to be more flexible
a = 0.2;
c0 = 299792458;
f = 2*10^9; % 2 ghz
k0 = (f/c0)^2 * 4*pi^2;

% only needed for small b
k_z10 = sqrt(k0^2-(pi/a)^2); 
%Einc = e10*
%gamma = 1j*k_z10;


% Quadrature rule
q2u = [[6.666666666666667e-01, ...
        1.666666666666667e-01]; ...
       [1.666666666666667e-01, ...
        6.666666666666667e-01]; ...
       [1.666666666666667e-01, ...
        1.666666666666667e-01]]';
q2w = [1.666666666666667e-02; ...
       1.666666666666667e-02; ...
       1.666666666666667e-02]';

% H(grad) basis functions
% Node basis functions for the reference element
up{1} = 1 - q2u(1,:) - q2u(2,:);
up{2} = q2u(1,:);
up{3} = q2u(2,:);

% Gradient of H(grad) basis functions, used in H(div), H(curl)
ug{1} = [-1 -1]';
ug{2} = [+1 0]';
ug{3} = [0 +1]';

% H(curl) basis function N
uin{1} = [ug{2}*up{1} - ug{1}*up{2}; 0,0,0];
uin{2} = [ug{3}*up{2} - ug{2}*up{3}; 0,0,0];
uin{3} = [ug{1}*up{3} - ug{3}*up{1}; 0,0,0];

%S, n x N
usn{1} = [-up{1}*ug{2}(2)+up{2}*ug{1}(2); up{1}*ug{2}(1)-up{2}*ug{1}(1); 0,0,0];
usn{2} = [-up{2}*ug{3}(2)+up{2}*ug{3}(2); up{2}*ug{3}(1)-up{2}*ug{3}(1); 0,0,0];
usn{3} = [-up{3}*ug{1}(2)+up{1}*ug{3}(2); up{3}*ug{1}(1)-up{1}*ug{3}(1); 0,0,0];

usnn{1} = [-up{2}*ug{1}(1)+up{1}*ug{2}(1); -up{2}*ug{1}(2)+up{1}*ug{2}(2); 0,0,0];
usnn{2} = [-up{2}*ug{3}(1)+up{2}*ug{3}(1); -up{2}*ug{3}(2)+up{2}*ug{3}(2); 0,0,0];
usnn{3} = [-up{3}*ug{1}(1)+up{1}*ug{3}(1); -up{3}*ug{1}(2)+up{1}*ug{3}(2); 0,0,0];



% Physical coordinates
% Maps from refference element to physical element
q2x = zeros(3,length(q2w)); % allocate memorry
for iIdx = 1:3 % 3 edge points i.e triangle
    q2x = q2x + xyz(:,iIdx)*up{iIdx};
end

% Jacobian
jac = zeros(3);
for iIdx = 1:3
    jac = jac ...
        + [xyz(1,iIdx)*ug{iIdx}', 0; ...
           xyz(2,iIdx)*ug{iIdx}', 0; ...
           xyz(3,iIdx)*ug{iIdx}', 1]';
end

% Mappings
det_jac = det(jac);
map_ccs = inv(jac);      % mapping for curl-conforming space
%map_dcs = jac'/det_jac;  % mapping for div-conforming space



% Jacobian for U_inc
jacU = zeros(3);
for iIdx = 1:3
    jac = jac ...
        + [xyz(1,iIdx)*ug{iIdx}', 0; ...
           xyz(2,iIdx)*ug{iIdx}', 0; ...
           xyz(3,iIdx)*ug{iIdx}', 1]';
end

% Mappings
det_jac = det(jac);
map_ccs = inv(jac);




for iIdx = 1:3
    gin{iIdx} = map_ccs*uin{iIdx};
    gsn{iIdx} = map_ccs*usn{iIdx};
    gsnn{iIdx} = map_ccs*usnn{iIdx};
end

% B_{ij} ElMtx [j^{-1}] n x < N
for iIdx = 1:3
    for jIdx = 1:3
        %maTmp = ones(size(q2w));
        ipTmp = sum(gsn{iIdx}.* gsn{jIdx});
        BElMtx_EE(iIdx,jIdx) = gamma * ipTmp * q2w' * det_jac;
    end
end

%b_{ij} ElMtx n x n x [j^{-1}]N
e10xy = sin()
for iIdx = 1:3
    ipTmp = sum(gsnn{iIdx});
    
    bElMtx_EE(iIdx,jIdx) = gamma * ipTmp * q2w' * det_jac;
end

