function [T_db] = ...
    S_parameters(E,ed2no_port1,ed2no_port2,fac2no_port1,fac2no_port2,no2xyz,a) % need to change from hardcode to adaptive

% ed2no_port1
% fac2no_port1

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

direction = 1;
n = repmat([0,0,1],3,1)'*direction;

%S, n x N
usn{1} = cross(n,uin{1});
usn{2} = cross(n,uin{2});
usn{3} = cross(n,uin{3});

%n x S (n x (n x S)) 
usnn{1} = cross(n,usn{1});
usnn{2} = cross(n,usn{2});
usnn{3} = cross(n,usn{3});

res = 0;

figure(15), clf, hold on
for no = fac2no_port2
    
    xyz = no2xyz(:,no);
    
    noTmp = zeros(2,3); %zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
    noTmp(:) = no([[1,2];[2,3];[3,1]]');%el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
    esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
    eiVtr = ElementDatabase_Get('edges', noTmp); % each edge id assosiated with the specific thetra (6 edges)

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
    map_ccs = inv(jac);   
    
    Eel = imag(E(eiVtr)).*esVtr';
    
    for iIdx = 1:3
        gsnn{iIdx} = map_ccs*usn{iIdx};
    end
    
    e10 = [0,0,0;sin( (pi* (q2x(1,:)+(a/2))) / a );0,0,0];
    
    
    ipTmp = zeros(3,3);
    
    %go throug each basis function
    plot([xyz(1,:), xyz(1,1)],[xyz(2,:),xyz(2,1)],'Color',[0.35, 0.35, 0.35])
    scale = 2;
    for basis = 1:3
        % go thrug each point with in the 3 sampled points in the triangle
        tmp = gsnn{basis}*Eel(iIdx)*det_jac;
        quiver(q2x(1,:),q2x(2,:),tmp(1,:)*scale,tmp(2,:)*scale,'AutoScale','off')
        ipTmp = ipTmp + tmp;
    end
    
    
    quiver(q2x(1,:),q2x(2,:),ipTmp(1,:)*scale,ipTmp(2,:)*scale,'k','AutoScale','off')
    
    %ipTmp = ipTmp/3
    
    %multiply each point with e10
    for point = 1:3
       magTemp(point) = ipTmp(:,point)'*e10(:,point);
    end
    
    res = res + magTemp* q2w' %* det_jac;
    
end
    
%constand before integral, for both R and T
nrom = @(z) 2*exp(-1j*k_z10*z);
z_val = no2xyz(:,fac2no_port2(:));
z_val = mean(z_val(3,:));
T = norm(z_val)*res;
T_db = 10*log10(T);

