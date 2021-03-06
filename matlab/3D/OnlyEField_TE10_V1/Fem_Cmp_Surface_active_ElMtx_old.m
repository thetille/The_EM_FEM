% --------------------------------------------------------------
% Compute element matrices for the port surfaces (trangles) by means of
% numerical integration on the reference element
% --------------------------------------------------------------
function [bElMtx_EE] = ...
    Fem_Cmp_Surface_ElMtx(xyz,k_z10,a,E0) % need to change from hardcode to adaptive
% Argument:
%   xyz = the coordinates of the nodes of the element
%   ma2er = material to permittivity
%   ma2si = material to conductivity
% Return:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global normals;
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

direction = -1;
n = repmat([0,0,1],3,1)'*direction;
if normals
    figure(3);
    scale = 0.05;
    quiver3(xyz(1,:),xyz(2,:),xyz(3,:),n(1,:)*scale,n(2,:)*scale,n(3,:)*scale,'Autoscale', 'off')
end


% usn{1} = [-up{1}*ug{2}(2)+up{2}*ug{1}(2); up{1}*ug{2}(1)-up{2}*ug{1}(1); 0,0,0];
% usn{2} = [-up{2}*ug{3}(2)+up{3}*ug{2}(2); up{2}*ug{3}(1)-up{3}*ug{2}(1); 0,0,0];
% usn{3} = [-up{3}*ug{1}(2)+up{1}*ug{3}(2); up{3}*ug{1}(1)-up{1}*ug{3}(1); 0,0,0];


% usnn{1} = [-up{1}*ug{2}(1)+up{2}*ug{1}(1); -up{1}*ug{2}(2)+up{2}*ug{1}(2); 0,0,0];
% usnn{2} = [-up{2}*ug{3}(1)+up{3}*ug{2}(1); -up{2}*ug{3}(2)+up{3}*ug{2}(2); 0,0,0];
% usnn{3} = [-up{3}*ug{1}(1)+up{1}*ug{3}(1); -up{3}*ug{1}(2)+up{1}*ug{3}(2); 0,0,0];



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
det_jac = det(jac);
map_ccs = inv(jac);       % mapping for curl-conforming space
%map_dcs = jac'/det_jac;  % mapping for div-conforming space



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


%b_{ij} ElMtx n x n x [j^{-1}]N * U_inc
% %Jacobian for U_inc
% jacU = zeros(3);
% jacU = [xyz(1,3)-xyz(1,3), xyz(1,2)-xyz(1,1); ...
%         xyz(2,3)-xyz(2,1), xyz(2,2)-xyz(2,1)];
% % Mappings
% det_jacU = det(jacU);
% 
% figure(10), hold on
% scatter(xyz(1,:),xyz(2,:),'g');
% scatter(q2x(1,:),q2x(2,:),'r','x');

%xyz = xyz+(a/2); % (the waveguide needs to start at 0 to a
%Area = 1/2*abs(xyz(1,1)*(xyz(2,2)-xyz(2,3))+ xyz(1,2)*(xyz(2,3)-xyz(2,1)) + xyz(1,3)*(xyz(2,1)-xyz(2,2)));

offset = (a/2);
Z = mean(xyz(3,:));

Uinc = -2j*k_z10*E0*[0,0,0;...
                     sin((pi*(q2x(1,:)+offset)) / a );...
                     0,0,0]...
                     .*exp(-1j*k_z10*Z); %needs z in exponent

% quiver(q2x(1,:),q2x(2,:),Uinc(1,:),imag(Uinc(2,:)*scale),'Autoscale', 'off')
%Uinc_Int = Uinc_Int + Uinc(2,:)*q2w'*det_jac;

for iIdx = 1:3
    ipTmp = sum(usnn{iIdx} .* Uinc);
    bElMtx_EE(iIdx) = (ipTmp * q2w')*det_jac;
end

% 
% 
% % ipTmp = [usnn{1}(:),usnn{2}(:),usnn{3}(:)]' * Uinc(:);
% % bElMtx_EE(iIdx) = (ipTmp' * q2w')*det_jac ;
% 
% scale = 0.001;
% %  figure(20), hold on
% %  %plot triangl
% %  plot([xyz(1,:) xyz(1,1)],[xyz(2,:) xyz(2,1)])
% %  %plot basis functions
% % for iIdx = 1:3
% %     quiver(q2x(1,:),q2x(2,:),usnn{iIdx}(1,:)*imag(bElMtx_EE(iIdx))*scale,usnn{iIdx}(2,:)*imag(bElMtx_EE(iIdx))*scale,'Autoscale', 'off')
% % end
% 
% figure(21), hold on
% %plot triangl
% plot([xyz(1,:) xyz(1,1)],[xyz(2,:) xyz(2,1)])
% %plot basis functions
% x_dir = 0;
% y_dir = 0;
% 
% for iIdx = 1:3
%     x_dir = x_dir + gin{iIdx}(1,:)*imag(bElMtx_EE(iIdx));
%     y_dir = y_dir + gin{iIdx}(2,:)*imag(bElMtx_EE(iIdx));
% end
% 
% quiver(q2x(1,:),q2x(2,:),x_dir*scale,y_dir*scale,'Autoscale', 'off')
% 
% X = [X, q2x(1,:)];
% Y = [Y, q2x(2,:)];
% Z = [Z, y_dir.^2+x_dir.^2];
% 
% figure(23)
% 
% scale = det_jac*2000;
% 
% hold on
% plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
% daspect([1,1,1])
% start = q2x;
% % scale = 0.0001;
% for iIdx = 1:3
%     quiver(start(1,:),start(2,:), gin{iIdx}(1,:)*imag(bElMtx_EE(iIdx))*scale, gin{iIdx}(2,:)*imag(bElMtx_EE(iIdx))*scale,'Autoscale','off')
%     start(1,:) = start(1,:) + gin{iIdx}(1,:)*imag(bElMtx_EE(iIdx))*scale;
%     start(2,:) = start(2,:) + gin{iIdx}(2,:)*imag(bElMtx_EE(iIdx))*scale;
% end
% %scale = 0.1;
% 
% quiver(q2x(1,:),q2x(2,:),x_dir*scale,y_dir*scale,'Autoscale', 'off')
% quiver(q2x(1,:),q2x(2,:),Uinc(1,:),imag(Uinc(2,:)*det_jac*scale),'k','Autoscale', 'off')
% integralx = integralx + x_dir*q2w'*det_jac;
% integraly = integraly + y_dir*q2w'*det_jac;
% 
% % figure(22)
% % pgon = polyshape([xyz(1,:) xyz(1,1)],[xyz(2,:) xyz(2,1)]);
% % 
% % plot(pgon,'FaceColor','red','FaceAlpha',0.1)

