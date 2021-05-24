% --------------------------------------------------------------
% Compute element matrices for the port surfaces (trangles) by means of
% numerical integration on the reference element
% --------------------------------------------------------------
function [bElMtx_EE] = ...
    Fem_Cmp_Surface_ElMtx(xyz,k_z10,a,E0,direction) % need to change from hardcode to adaptive
% Argument:
%   xyz = the coordinates of the nodes of the element
%   ma2er = material to permittivity
%   ma2si = material to conductivity
% Return:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global normals;
global temp1;
global temp2;
normals = 0;

% Quadrature rule
q2u = [[6.666666666666667e-01, ...
        1.666666666666667e-01]; ...
       [1.666666666666667e-01, ...
        6.666666666666667e-01]; ...
       [1.666666666666667e-01, ...
        1.666666666666667e-01]]';
q2w = [1.666666666666667e-01; ...
       1.666666666666667e-01; ...
       1.666666666666667e-01]';

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

%direction = direction*-1;
n = repmat([0,0,1],3,1)'*direction;

if normals
    figure(3);
    scale = 0.05;
    quiver3(xyz(1,:),xyz(2,:),xyz(3,:),n(1,:)*scale,n(2,:)*scale,n(3,:)*scale,'Autoscale', 'off')
end

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
           xyz(3,iIdx)*ug{iIdx}', 0]';
end
jac(3,3) = 1;
% Mappings
map_ccs = inv(jac); % mapping for curl-conforming space

for iIdx = 1:3
    gin{iIdx} = map_ccs*uin{iIdx};
end
%S, n x N
usn{1} = cross(n,gin{1});
usn{2} = cross(n,gin{2});
usn{3} = cross(n,gin{3});


%n x S (n x (n x S)) 
usnn{1} = cross(n,usn{1});
usnn{2} = cross(n,usn{2});
usnn{3} = cross(n,usn{3});

zval = mean(xyz(3,:));

% fprintf("-2j*k_z10*E0: %f\n",-2j*k_z10*E0)
% fprintexf("exp(-1j*k_z10*zval): %f\n",exp(-1j*k_z10*zval))

Uinc = -2j*k_z10*E0*[0,0,0;sin((pi*(q2x(1,:)+(a/2))) / a );0,0,0].*exp(-1j*k_z10*zval);%*direction; %needs z in exponent

plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
quiver(q2x(1,:),q2x(2,:),abs(Uinc(1,:)),abs(Uinc(2,:)))


% clf
% quiver(q2x(1,:),q2x(2,:),Uinc(1,:),Uinc(2,:))
% 
% temp1 = temp1 + Uinc*q2w'*det(jac);
% fprintf('Uinc: %f + %fi\n',real(temp1(2)),imag(temp1(2)))

%bElMtx_EE = ([usnn{1}(:), usnn{2}(:), usnn{3}(:)]\(Uinc(:)*det(jac)))*direction;


for iIdx = 1:3

    ipTmp = sum(usnn{iIdx} .* Uinc);
    bElMtx_EE(iIdx) = (ipTmp * q2w')*det(jac);
end



% uincsim = 0;
% for iIdx = 1:3
%     uincsim = uincsim + usnn{iIdx}*bElMtx_EE(iIdx);
% end
% hold on
% quiver(q2x(1,:),q2x(2,:),uincsim(1,:),uincsim(2,:))
% 
% temp2 = temp2+ uincsim*q2w';
% 
% 
% fprintf('this: %f + %fi\n',real(temp2(2)),imag(temp2(2)))