clear
close all
a = 10;
%xyz = [1,2,3.5;3,1,2;0,0,0];
xyz = [rand(2,3)*a; 0,0,0];
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
scale = 0.2;
figure(1), clf
subplot(3,1,1), hold on
plot([0,0,1,0]',[0,1,0,0]')
scatter(q2u(1,:),q2u(2,:))
quiver(q2u(1,:),q2u(2,:),uin{1}(1,:)*scale,uin{1}(2,:)*scale,'Autoscale','off')

subplot(3,1,2), hold on
plot([0,0,1,0]',[0,1,0,0]')
scatter(q2u(1,:),q2u(2,:))
quiver(q2u(1,:),q2u(2,:),uin{2}(1,:)*scale,uin{2}(2,:)*scale,'Autoscale','off')

subplot(3,1,3), hold on
plot([0,0,1,0]',[0,1,0,0]')
scatter(q2u(1,:),q2u(2,:))
quiver(q2u(1,:),q2u(2,:),uin{3}(1,:)*scale,uin{3}(2,:)*scale,'Autoscale','off')


figure(2), clf, hold on
plot([0,0,a,a,0],[0,a,a,0,0])
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')


for iIdx = 1:3
    gsnn{iIdx} = map_ccs*uin{iIdx};
end

scale = det_jac*0.2;
figure(3)
subplot(3,1,1), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),gsnn{1}(1,:)*scale,gsnn{1}(2,:)*scale,'Autoscale','off')
daspect([1 1 1])

subplot(3,1,2), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),gsnn{2}(1,:)*scale,gsnn{2}(2,:)*scale,'Autoscale','off')
daspect([1 1 1])

subplot(3,1,3),hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),gsnn{3}(1,:)*scale,gsnn{3}(2,:)*scale,'Autoscale','off')
daspect([1 1 1])

n = repmat([0,0,1],3,1)';
%S, n x N
usn{1} = cross(n,gsnn{1});
usn{2} = cross(n,gsnn{2});
usn{3} = cross(n,gsnn{3});

col = [0.75,0.75,0.75];

figure(4)
subplot(3,1,1), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),gsnn{1}(1,:)*scale,gsnn{1}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usn{1}(1,:)*scale,usn{1}(2,:)*scale,'Autoscale','off')
daspect([1 1 1])

subplot(3,1,2), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),gsnn{2}(1,:)*scale,gsnn{2}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usn{2}(1,:)*scale,usn{2}(2,:)*scale,'Autoscale','off')
daspect([1 1 1])

subplot(3,1,3),hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),gsnn{3}(1,:)*scale,gsnn{3}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usn{3}(1,:)*scale,usn{3}(2,:)*scale,'Autoscale','off')
daspect([1 1 1])

%n x S (n x (n x S)) 
usnn{1} = cross(n,usn{1});
usnn{2} = cross(n,usn{2});
usnn{3} = cross(n,usn{3});

figure(5)
subplot(3,1,1), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),usn{1}(1,:)*scale,usn{1}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usnn{1}(1,:)*scale,usnn{1}(2,:)*scale,'Autoscale','off')
daspect([1 1 1])

subplot(3,1,2), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),usn{2}(1,:)*scale,usn{2}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usnn{2}(1,:)*scale,usnn{2}(2,:)*scale,'Autoscale','off')
daspect([1 1 1])

subplot(3,1,3),hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),usn{3}(1,:)*scale,usn{3}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usnn{3}(1,:)*scale,usnn{3}(2,:)*scale,'Autoscale','off')
daspect([1 1 1])


%%%%%% old version %%%%%%%

% usn_old{1} = cross(n,uin{1});
% usn_old{2} = cross(n,uin{2});
% usn_old{3} = cross(n,uin{3});

% for iIdx = 1:3
%     gsn_old{iIdx} = map_ccs*usn_old{iIdx};
% end
% 
% figure(6)
% subplot(3,1,1), hold on
% plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
% scatter(q2x(1,:),q2x(2,:),'x')
% quiver(q2x(1,:),q2x(2,:),uin{1}(1,:)*scale,uin{1}(2,:)*scale,'Autoscale','off','color',col)
% quiver(q2x(1,:),q2x(2,:),gsn_old{1}(1,:)*scale,gsn_old{1}(2,:)*scale,'Autoscale','off')
% daspect([1 1 1])
% 
% subplot(3,1,2), hold on
% plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
% scatter(q2x(1,:),q2x(2,:),'x')
% quiver(q2x(1,:),q2x(2,:),uin{2}(1,:)*scale,uin{2}(2,:)*scale,'Autoscale','off','color',col)
% quiver(q2x(1,:),q2x(2,:),gsn_old{2}(1,:)*scale,gsn_old{2}(2,:)*scale,'Autoscale','off')
% daspect([1 1 1])
% 
% subplot(3,1,3),hold on
% plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
% scatter(q2x(1,:),q2x(2,:),'x')
% quiver(q2x(1,:),q2x(2,:),uin{3}(1,:)*scale,uin{3}(2,:)*scale,'Autoscale','off','color',col)
% quiver(q2x(1,:),q2x(2,:),gsn_old{3}(1,:)*scale,gsn_old{3}(2,:)*scale,'Autoscale','off')
% daspect([1 1 1])


figure(7)

Uinc = [0,0,0;sin((pi*(q2x(1,:)))/a);0,0,0];
UincScaled = Uinc*det_jac;

for iIdx = 1:3
    ipTmp = sum(usnn{iIdx} .* UincScaled);%*det_jac;
    bElMtx_EE(iIdx) = (ipTmp * q2w') ;
end
hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),usnn{3}(1,:)*scale,usnn{3}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),Uinc(1,:)*scale,Uinc(2,:)*scale,'Autoscale','off')

%quiver(q2x(1,:),q2x(2,:),usnn{3}(1,:).*ipTmp{3}*scale,usnn{3}(2,:).*ipTmp{3}*scale,'Autoscale','off')
daspect([1 1 1])


figure(8), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
daspect([1,1,1])
start = q2x;
for iIdx = 1:3
    quiver(start(1,:),start(2,:),usnn{iIdx}(1,:)*bElMtx_EE(iIdx)*scale,usnn{iIdx}(2,:).*bElMtx_EE(iIdx)*scale,'Autoscale','off')
    start(1,:) = start(1,:) + usnn{iIdx}(1,:).*bElMtx_EE(iIdx)*scale;
    start(2,:) = start(2,:) + usnn{iIdx}(2,:).*bElMtx_EE(iIdx)*scale;
end
% quiver(q2x(1,:),q2x(2,:),usnn{1}(1,:).*ipTmp{1}*scale,usnn{1}(2,:).*ipTmp{1}*scale,'Autoscale','off')
% quiver(q2x(1,:),q2x(2,:),usnn{2}(1,:).*ipTmp{2}*scale,usnn{2}(2,:).*ipTmp{2}*scale,'Autoscale','off')
% quiver(q2x(1,:),q2x(2,:),usnn{3}(1,:).*ipTmp{3}*scale,usnn{3}(2,:).*ipTmp{3}*scale,'Autoscale','off')


%figure(9), hold on
%plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
quiver(q2x(1,:),q2x(2,:),UincScaled(1,:)*scale*(-1),UincScaled(2,:)*scale*(-1),'k','Autoscale','off')


figure(9), hold on

et = [usnn{1}(:), usnn{2}(:), usnn{3}(:)]\UincScaled(:);


plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
daspect([1,1,1])
start = q2x;
for iIdx = 1:3
    quiver(start(1,:),start(2,:),usnn{iIdx}(1,:)*et(iIdx)*scale,usnn{iIdx}(2,:)*et(iIdx)*scale,'Autoscale','off')
    start(1,:) = start(1,:) + usnn{iIdx}(1,:)*et(iIdx)*scale;
    start(2,:) = start(2,:) + usnn{iIdx}(2,:)*et(iIdx)*scale;
end
% quiver(q2x(1,:),q2x(2,:),usnn{1}(1,:).*ipTmp{1}*scale,usnn{1}(2,:).*ipTmp{1}*scale,'Autoscale','off')
% quiver(q2x(1,:),q2x(2,:),usnn{2}(1,:).*ipTmp{2}*scale,usnn{2}(2,:).*ipTmp{2}*scale,'Autoscale','off')
% quiver(q2x(1,:),q2x(2,:),usnn{3}(1,:).*ipTmp{3}*scale,usnn{3}(2,:).*ipTmp{3}*scale,'Autoscale','off')


%figure(9), hold on
%plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
quiver(q2x(1,:),q2x(2,:),UincScaled(1,:)*scale,UincScaled(2,:)*scale,'k','Autoscale','off')


%b_{ij} ElMtx n x n x [j^{-1}]N * U_inc
% %Jacobian for U_inc
% jacU = zeros(3);
% jacU = [xyz(1,3)-xyz(1,3), xyz(1,2)-xyz(1,1); ...
%         xyz(2,3)-xyz(2,1), xyz(2,2)-xyz(2,1)];
% % Mappings
% det_jacU = det(jacU);

% figure(10), hold on
% scatter(xyz(1,:),xyz(2,:),'g');
% scatter(q2x(1,:),q2x(2,:),'r','x');

%xyz = xyz+(a/2); % (the waveguide needs to start at 0 to a
%Area = 1/2*abs(xyz(1,1)*(xyz(2,2)-xyz(2,3))+ xyz(1,2)*(xyz(2,3)-xyz(2,1)) + xyz(1,3)*(xyz(2,1)-xyz(2,2)));

%Uinc = -2j*k_z10*E0*[0,0,0;sin((pi*(q2x(1,:)+(a/2))) / a );0,0,0].*exp(-1j*k_z10*0); %needs z in exponent
% % scale = 0.001;
% % quiver(q2x(1,:),q2x(2,:),Uinc(1,:),imag(Uinc(2,:)*scale),'Autoscale', 'off')
% 
% for iIdx = 1:3
%     ipTmp = sum(usnn{iIdx} .* Uinc);
%     bElMtx_EE(iIdx) = (ipTmp * q2w') *det_jac;
% end
