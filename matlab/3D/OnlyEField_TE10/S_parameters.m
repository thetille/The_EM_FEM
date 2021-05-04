function [S_par] = ...
    S_parameters(E,fac2no_port1,fac2no_port2,no2xyz,a,k_z10) % need to change from hardcode to adaptive

% ed2no_port1
% fac2no_port1

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

% direction = 1;
% n = repmat([0,0,1],3,1)'*direction;
% 
% %S, n x N
% usn{1} = cross(n,uin{1});
% usn{2} = cross(n,uin{2});
% usn{3} = cross(n,uin{3});
% 
% %n x S (n x (n x S)) 
% usnn{1} = cross(n,usn{1});
% usnn{2} = cross(n,usn{2});
% usnn{3} = cross(n,usn{3});


port_list{1} = fac2no_port1;
port_list{2} = fac2no_port2;

for port = 1:length(port_list)
    res = 0;
    for no = port_list{port}

        xyz = no2xyz(:,no);
        noTmp = zeros(2,3); %zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
        noTmp(:) = no([[1,2];[2,3];[3,1]]');%el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
        esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
        eiVtr = ElementDatabase_Get('edges', noTmp); % each edge id assosiated with the specific triangle (3 edges)
        Eel = E(eiVtr).*esVtr';

        % Physical coordinates
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

        gin = cell(3,1);
        for iIdx = 1:3
            gin{iIdx} = jac\uin{iIdx};  % sam as map_ccs = inv(jac); , uni{iIdx} * map_css
        end

        e10 = [0,0,0;sin( (pi* (q2x(1,:)+(a/2))) / a );0,0,0];
        %go throug each basis function
        ipTmp = zeros(3,3);
        for basis = 1:3
            % go thrug each point with in the 3 sampled points in the triangle
            tmp = gin{basis}*Eel(basis);
            ipTmp = ipTmp + tmp;

        end
        


        %multiply each point with e10
        magTemp = zeros(1,3);
        for point = 1:3
           magTemp(point) = ipTmp(:,point)'*e10(:,point)*0.5; % the 0.5 is not part of the equation and it would be best if it was not there.
        end
        res = res + magTemp * q2w' * det_jac;
         
%       scale = 10;
%       figure(15+port), clf, hold on,
%       plot(q2x(1,:),q2x(2,:),'x')
%       plot([xyz(1,:), xyz(1,1)],[xyz(2,:),xyz(2,1)],'Color',[0.35, 0.35, 0.35])
%       quiver(q2x(1,:),q2x(2,:),e10(1,:)*scale*det_jac,e10(2,:)*scale*det_jac,'b','AutoScale','off')
%       quiver(q2x(1,:),q2x(2,:),abs(ipTmp(1,:)*det_jac)*scale,abs(ipTmp(2,:)*det_jac)*scale,'k','AutoScale','off')
    end
    %constand before integral, for both R and T
    TRConstant = @(z) (2/(a*a));
    z_val = no2xyz(:,fac2no_port2(:));
    z_val = mean(z_val(3,:));
    S_par(port) = (TRConstant(z_val)*res)^2;
end
%S_par(:,1) = S_par(:,1)*(-1)-1;

