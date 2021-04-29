function [T_db] = ...
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

res = 0;
res2 = 0;
plot_debug = 1;
test_res = 0;
test_res2 = 0;
test_res3 = 0;
if plot_debug, figure(15), clf, hold on, end
for no = fac2no_port2
    
    xyz = no2xyz(:,no);
    
    noTmp = zeros(2,3); %zeros(size(ed2noLoc)); % temp var of size of the amount of initial base lines (edges) 
    noTmp(:) = no([[1,2];[2,3];[3,1]]');%el2no(ed2noLoc(:),elIdx); % temp var with initial base edge nodes (one edge has 2 nodes ie 6x2)
    esVtr = sign(noTmp(2,:)-noTmp(1,:)); % each edge gest an assigned direction based on node id (this is random)is only used as initial value
    eiVtr = ElementDatabase_Get('edges', noTmp); % each edge id assosiated with the specific triangle (3 edges)

    % Physical coordinates
    % Maps from refference element to physical element
    q2x = zeros(3,length(q2w)); % allocate memorry
    for iIdx = 1:3 % 3 edge points i.e triangle
        q2x = q2x + xyz(:,iIdx)*up{iIdx};
    end
    
    if plot_debug, plot(q2x(1,:),q2x(2,:),'x'), end
    
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
    
    Eel = E(eiVtr).*esVtr';
    
    for iIdx = 1:3
        gin{iIdx} = map_ccs*uin{iIdx};
    end
    
    % Evaluation of element matrix: phi_i phi_j
    iElMtx_NN = zeros(3);
    for iIdx = 1:3
        for jIdx = 1:3
            bsiTmp = up{iIdx};
            bsjTmp = up{jIdx};
            iElMtx_NN(iIdx,jIdx) = (bsiTmp.*bsjTmp) * q2w';
        end
    end

    % Evaluation of element matrix: phi_i (x * Nj)
    iElMtx_NE.x = zeros(3,3);
    iElMtx_NE.y = zeros(3,3);
    iElMtx_NE.z = zeros(3,3);
    for iIdx = 1:3
        for jIdx = 1:3
            bsiTmp = up{iIdx};
            bxjTmp = gin{jIdx}(1,:);
            byjTmp = gin{jIdx}(2,:);
            bzjTmp = gin{jIdx}(3,:);
            iElMtx_NE.x(iIdx,jIdx) = (bsiTmp.*bxjTmp) * q2w' * det_jac;
            iElMtx_NE.y(iIdx,jIdx) = (bsiTmp.*byjTmp) * q2w' * det_jac;
            iElMtx_NE.z(iIdx,jIdx) = (bsiTmp.*bzjTmp) * q2w' * det_jac;
        end
    end
    
    isMtx = diag(1./diag(iElMtx_NN));
% 
    Eelxc = isMtx*iElMtx_NE.x*Eel;
    Eelyc = isMtx*iElMtx_NE.y*Eel;
    Eelzc = isMtx*iElMtx_NE.z*Eel;

%     Eelxc = iElMtx_NE.x*Eel;
%     Eelyc = iElMtx_NE.y*Eel;
%     %Eelzc = iElMtx_NE.z*Eel;
    
    if plot_debug, plot([xyz(1,:), xyz(1,1)],[xyz(2,:),xyz(2,1)],'Color',[0.35, 0.35, 0.35]), end
    
    %scale = 400;
    %if plot_debug, quiver(q2x(1,:),q2x(2,:),abs(Eelxc)'*scale,abs(Eelyc)'*scale,'k','AutoScale','off'), end
    
    e10 = [0,0,0;sin( (pi* (q2x(1,:)+(a/2))) / a );0,0,0];
%    if plot_debug, quiver(q2x(1,:),q2x(2,:),e10(1,:)*scale,e10(2,:)*scale,'k','AutoScale','off'), end
    scale = 10;
    if plot_debug, quiver(q2x(1,:),q2x(2,:),e10(1,:)*scale*det_jac,e10(2,:)*scale*det_jac,'b','AutoScale','off'), end
%     
%     
     %Eelxc = Eelxc.*(e10(1,:)');
     %Eelyc = Eelyc.*(e10(2,:)');
     
     test_res = test_res + e10(2,:) *det_jac * q2w' ;
     %test_res2 = test_res2 + abs(Eelyc)' * q2w' ;
     
%     Area = 1/2*abs(xyz(1,1)*(xyz(2,2)-xyz(2,3))+ xyz(1,2)*(xyz(2,3)-xyz(2,1)) + xyz(1,3)*(xyz(2,1)-xyz(2,2)));
%      disp(Area)
%      disp(det_jac)
%      disp(det_jac/Area)
%     
     %res = res + Eelyc' * q2w'*det_jac;
     %res2 = res2 + sum(Eelyc)*det_jac;

     
    ipTmp = zeros(3,3);
    %go throug each basis function
    scale = 10;
    for basis = 1:3
        % go thrug each point with in the 3 sampled points in the triangle
        tmp = gin{basis}*Eel(basis);
        %quiver(q2x(1,:),q2x(2,:),tmp(1,:)*scale,tmp(2,:)*scale,'AutoScale','off')
        ipTmp = ipTmp + tmp;
    end
    quiver(q2x(1,:),q2x(2,:),abs(ipTmp(1,:)*det_jac)*scale,abs(ipTmp(2,:)*det_jac)*scale,'k','AutoScale','off')
    
    %ipTmp = ipTmp/3
    
%     multiply each point with e10
%     for point = 1:3
%        magTemp(point) = ipTmp(:,point)'*e10(:,point);
%     end
%     
     test_res3 = test_res3 + abs(ipTmp(2,:))* q2w' * det_jac;
    
end
fprintf("test_res: %f\n",abs(test_res))
fprintf("test_res2: %f\n",abs(test_res2))
fprintf("test_res3: %f\n",abs(test_res3))
%constand before integral, for both R and T
TRConstant = @(z) (2*exp(-1j*k_z10*z))/(a*a);
z_val = no2xyz(:,fac2no_port2(:));
z_val = mean(z_val(3,:));
T = abs(TRConstant(z_val)*res);
T2 = abs(TRConstant(z_val)*res2);
%T3 = abs(TRConstant(z_val)*res3);

T12 = T^2;
T22 = T2^2;
%T32 = T3^2;

fprintf("T: %f \t T^2: %f\n",T,T12)
fprintf("T2: %f \t T2^2: %f\n",T2,T22)

xlim([-0.1,0.1])
ylim([-0.1,0.1])
pbaspect([1,1,1])
%daspect([1 1 1])
disp('hej')
%T_db = 10*log10(T2)


