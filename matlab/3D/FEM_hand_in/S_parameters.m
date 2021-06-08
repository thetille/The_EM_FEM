function [S_par] = ...
    S_parameters(E,port_fac2no_list,no2xyz,a,b,k_z10,E0) % need to change from hardcode to adaptive

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



for port = 1:length(port_fac2no_list)
    res = 0;
    int1 = 0;
    
    for no = port_fac2no_list{port}

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
        map_ccs = inv(jac);

        gin = cell(3,1);
        for iIdx = 1:3
            gin{iIdx} = map_ccs*uin{iIdx};  % sam as map_ccs = inv(jac); , uni{iIdx} * map_css
        end

        e10 = [0,0,0;sin( (pi* (q2x(1,:)+(a/2))) / a );0,0,0];
        %go throug each basis function
        ipTmp = zeros(3,3);
        for basis = 1:3
            % go thrug each point with in the 3 sampled points in the triangle
            tmp = gin{basis}*Eel(basis);
            ipTmp = ipTmp + tmp;

        end
        int1 = int1 + ipTmp * q2w' * det_jac;

        %multiply each point with e10
        magTemp = zeros(3,3);
        for point = 1:3
           magTemp(:,point) = ipTmp(:,point).*e10(:,point);
        end
        
        res = res + magTemp * q2w' * det_jac;
    end



    %S_par(port) = int1(2,:);
    S_par(port) = res(2);
    port_nodes = port_fac2no_list{port};
    temp = no2xyz(:,port_nodes(:));
    z_val(port) = mean(temp(3,:));
end

TConstant = @(z) ( (2*exp(1j*k_z10*z)) / (a*b*E0) );
RConstant = @(z) ( (2*exp(-1j*k_z10*z)) / (a*b*E0) );
S_par(:,1) = RConstant(z_val(1))*S_par(:,1)+exp(-2j*k_z10*z_val(1));
S_par(:,2) = TConstant(z_val(2))*S_par(:,2);
%S_par(:,3) = TConstant(z_val(3))*S_par(:,3);
exp(2j*k_z10*z_val(1))