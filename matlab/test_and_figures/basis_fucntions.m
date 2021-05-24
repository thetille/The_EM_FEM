clear
close all
a = 10;
xyz = [1,2,3.5;3,1,2;0,0,0];
%xyz = [rand(2,3)*a; 0,0,0];
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
   
   
res = 0.099;
point = [0,0];
u = [];
v = [];
while point(1) < 1
    while point(2) < (1-point(1))
        u = [u, point(1)];
        v = [v, point(2)];
        point(2) = point(2) + res;
    end
    point(1) = point(1) + res;
    point(2) = 0;
end

% H(grad) basis functions
% Node basis functions for the reference element
up{1} = 1 - q2u(1,:) - q2u(2,:);
up{2} = q2u(1,:);
up{3} = q2u(2,:);

upL{1} = 1-u-v;
upL{2} = u;
upL{3} = v;

% Gradient of H(grad) basis functions, used in H(div), H(curl)
ug{1} = [-1 -1]';
ug{2} = [+1 0]';
ug{3} = [0 +1]';

ugL{1} = [-1*ones(length(u),1), -1*ones(length(v),1)]';
ugL{2} = [1*ones(length(u),1) , 0*ones(length(v),1)]';
ugL{3} = [0*ones(length(u),1) , 1*ones(length(v),1)]';

% H(curl) basis function N
uin{1} = [ug{2}*up{1} - ug{1}*up{2}; 0,0,0];
uin{2} = [ug{3}*up{2} - ug{2}*up{3}; 0,0,0];
uin{3} = [ug{1}*up{3} - ug{3}*up{1}; 0,0,0];

uinL{1} = [ugL{2}.*upL{1} - ugL{1}.*upL{2};zeros(1,66)];
uinL{2} = [ugL{3}.*upL{2} - ugL{2}.*upL{3};zeros(1,66)];
uinL{3} = [ugL{1}.*upL{3} - ugL{3}.*upL{1};zeros(1,66)];

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
map_ccs2d = map_ccs(1:2,1:2);
%map_dcs = jac'/det_jac;  % mapping for div-conforming space




%% the 
col = [0.75,0.75,0.75];
col2 = [0.66,0.18,0.08];
col3 = [0,0,0];

scale = 0.2;
figure(1), clf, hold on
%subplot(1,1,1), hold on
h1 = plot([0,0,1,0]',[0,1,0,0]');
h2 = scatter(q2u(1,:),q2u(2,:),'x');
quiver(q2u(1,:),q2u(2,:),uin{1}(1,:)*scale,uin{1}(2,:)*scale,'Autoscale','off','color',col3,'DisplayName','N for quadrature rule')
quiver(u,v,uinL{1}(1,:)*scale,uinL{1}(2,:)*scale,'Autoscale', 'off','color',col2,'DisplayName','N')
ylim([-0.1,1.1])
daspect([1,1,1])

set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
legend()



%% trasnformed version

for iIdx = 1:3
    gsnn{iIdx} = map_ccs*uin{iIdx};
end

for iIdx = 1:3
    gsnnL{iIdx} = map_ccs*uinL{iIdx};
end

q2xL = zeros(2,length(u)); % allocate memorry
for iIdx = 1:3 % 3 edge points i.e triangle
    q2xL = q2xL + xyz(1:2,iIdx)*upL{iIdx};
end
figure(2), clf, hold on
h1 = plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)]); %plot triangle
h2 = scatter(q2x(1,:),q2x(2,:),'x');
quiver(q2x(1,:),q2x(2,:),gsnn{1}(1,:)*scale,gsnn{1}(2,:)*scale,'Autoscale','off','color',[0,0,0],'DisplayName','N for quadrature rule')
quiver(q2xL(1,:),q2xL(2,:),gsnnL{1}(1,:)*scale,gsnnL{1}(2,:)*scale,'Autoscale', 'off','color',col2,'DisplayName','N')
xlim([min(xyz(1,:))-0.1,max(xyz(1,:))+0.1])
ylim([min(xyz(2,:))-0.1,max(xyz(2,:))+0.1])
legend()

set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
legend()

%%
figure(3), clf, hold on
plot([0,0,a,a,0],[0,a,a,0,0])
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')

%scale = det_jac*0.2;
figure(4)
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
nL = repmat([0,0,1],66,1)';
n2 = repmat([0,0,1],3,1)'*-1;
%S, n x N
usn{1} = cross(n,gsnn{1});
usn{2} = cross(n,gsnn{2});
usn{3} = cross(n,gsnn{3});

usnL{1} = cross(nL,gsnnL{1});
usnL{2} = cross(nL,gsnnL{2});
usnL{3} = cross(nL,gsnnL{3});

usn2{1} = cross(n2,gsnn{1});
usn2{2} = cross(n2,gsnn{2});
usn2{3} = cross(n2,gsnn{3});


figure(5), hold on
h1 = plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)]);
h2 = scatter(q2x(1,:),q2x(2,:),'x');
quiver(q2xL(1,:),q2xL(2,:),gsnnL{1}(1,:)*scale,gsnnL{1}(2,:)*scale,'Autoscale','off','color',col,'DisplayName','N')
quiver(q2xL(1,:),q2xL(2,:),usnL{1}(1,:)*scale,usnL{1}(2,:)*scale,'Autoscale','off','DisplayName','S','color',col2)
quiver(q2x(1,:),q2x(2,:),usn{1}(1,:)*scale,usn{1}(2,:)*scale,'Autoscale','off','DisplayName','S for quadrature rule','color',col3)
daspect([1 1 1])
xlim([min(xyz(1,:))-0.1,max(xyz(1,:))+0.1])
ylim([min(xyz(2,:))-0.1,max(xyz(2,:))+0.1])

set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
legend()

%%
figure(6)
sgtitle('S')
subplot(2,3,1), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),gsnn{1}(1,:)*scale,gsnn{1}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usn{1}(1,:)*scale,usn{1}(2,:)*scale,'Autoscale','off')
quiver(q2x(1,:),q2x(2,:),usn2{1}(1,:)*scale,usn2{1}(2,:)*scale,'Autoscale','off','color',col2)
daspect([1 1 1])

subplot(2,3,2), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),gsnn{2}(1,:)*scale,gsnn{2}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usn{2}(1,:)*scale,usn{2}(2,:)*scale,'Autoscale','off')
quiver(q2x(1,:),q2x(2,:),usn2{2}(1,:)*scale,usn2{2}(2,:)*scale,'Autoscale','off','color',col2)
daspect([1 1 1])

subplot(2,3,3),hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),gsnn{3}(1,:)*scale,gsnn{3}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usn{3}(1,:)*scale,usn{3}(2,:)*scale,'Autoscale','off')
quiver(q2x(1,:),q2x(2,:),usn2{3}(1,:)*scale,usn2{3}(2,:)*scale,'Autoscale','off','color',col2)
daspect([1 1 1])

%n x S (n x (n x S)) 
usnn{1} = cross(n,usn{1});
usnn{2} = cross(n,usn{2});
usnn{3} = cross(n,usn{3});

usnnL{1} = cross(nL,usnL{1});
usnnL{2} = cross(nL,usnL{2});
usnnL{3} = cross(nL,usnL{3});

usnn2{1} = cross(n2,usn2{1});
usnn2{2} = cross(n2,usn2{2});
usnn2{3} = cross(n2,usn2{3});

figure(7), hold on
h1 = plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)]);
h2 = scatter(q2x(1,:),q2x(2,:),'x');
quiver(q2xL(1,:),q2xL(2,:),usnL{1}(1,:)*scale,usnL{1}(2,:)*scale,'Autoscale','off','color',col,'DisplayName','S')
quiver(q2xL(1,:),q2xL(2,:),usnnL{1}(1,:)*scale,usnnL{1}(2,:)*scale,'Autoscale','off','DisplayName','n x S','color',col2)
quiver(q2x(1,:),q2x(2,:),usnn{1}(1,:)*scale,usnn{1}(2,:)*scale,'Autoscale','off','DisplayName','n x S for quadrature rule','color',col3)
daspect([1 1 1])
xlim([min(xyz(1,:))-0.1,max(xyz(1,:))+0.1])
ylim([min(xyz(2,:))-0.1,max(xyz(2,:))+0.1])

set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
legend()

%figure(5)
%sgtitle('n x S')
subplot(2,3,4), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),usn{1}(1,:)*scale,usn{1}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usnn{1}(1,:)*scale,usnn{1}(2,:)*scale,'Autoscale','off')
quiver(q2x(1,:),q2x(2,:),usnn2{1}(1,:)*scale,usnn2{1}(2,:)*scale,'--','Autoscale','off','color',col2)
daspect([1 1 1])

subplot(2,3,5), hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),usn{2}(1,:)*scale,usn{2}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usnn{2}(1,:)*scale,usnn{2}(2,:)*scale,'Autoscale','off')
quiver(q2x(1,:),q2x(2,:),usnn2{2}(1,:)*scale,usnn2{2}(2,:)*scale,'--','Autoscale','off','color',col2)
daspect([1 1 1])

subplot(2,3,6),hold on
plot([xyz(1,:),xyz(1,1)],[xyz(2,:),xyz(2,1)])
scatter(q2x(1,:),q2x(2,:),'x')
quiver(q2x(1,:),q2x(2,:),usn{3}(1,:)*scale,usn{3}(2,:)*scale,'Autoscale','off','color',col)
quiver(q2x(1,:),q2x(2,:),usnn{3}(1,:)*scale,usnn{3}(2,:)*scale,'Autoscale','off','color',col2)
quiver(q2x(1,:),q2x(2,:),usnn2{3}(1,:)*scale,usnn2{3}(2,:)*scale,'--','Autoscale','off','color',col3)
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


figure(8)

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
