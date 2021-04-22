% see the transformation from ref to real for the basis elements
clear
%generate v,u cordinates
% 
res = 0.1;
point = [0,0];

q2u = [[6.666666666666667e-01, ...
        1.666666666666667e-01]; ...
       [1.666666666666667e-01, ...
        6.666666666666667e-01]; ...
       [1.666666666666667e-01, ...
        1.666666666666667e-01]]';

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


figure(1), clf
scatter(u,v); hold on
scatter(q2u(1,:),q2u(2,:),'g')


figure(2), clf
phi1 = 1-u-v;
phi2 = u;
phi3 = v;

phi1In = 1-q2u(1,:)-q2u(2,:);
phi2In = q2u(1,:);
phi3In = q2u(2,:);

scatter3(u,v,phi1)
sgtitle("phi")

scale = 0.05;
figure(3), clf
gPhi1 = [-1*ones(length(u),1) , -1*ones(length(v),1)]';
gPhi2 = [1*ones(length(u),1) , 0*ones(length(v),1)]';
gPhi3 = [0*ones(length(u),1) , 1*ones(length(v),1)]';

gPhi1In = [-1*ones(length(q2u(1,:)),1) , -1*ones(length(q2u(2,:)),1)]';
gPhi2In = [ 1*ones(length(q2u(1,:)),1)  , 0*ones(length(q2u(2,:)),1)]';
gPhi3In = [ 0*ones(length(q2u(1,:)),1)  , 1*ones(length(q2u(2,:)),1)]';

subplot(3,1,1)
quiver(u,v, gPhi1(1,:)*scale , gPhi1(2,:)*scale , 'Autoscale', 'off' ); hold on
quiver(q2u(1,:),q2u(2,:),gPhi1In(1,:)*scale,gPhi1In(2,:)*scale, 'Autoscale', 'off' ,'color','g')
subplot(3,1,2)
quiver(u,v,gPhi2(1,:)*scale,gPhi2(2,:)*scale, 'Autoscale', 'off' ), hold on
quiver(q2u(1,:),q2u(2,:),gPhi2In(1,:)*scale,gPhi2In(2,:)*scale, 'Autoscale', 'off' ,'color','g')
subplot(3,1,3)
quiver(u,v,gPhi3(1,:)*scale,gPhi3(2,:)*scale, 'Autoscale', 'off' ), hold on
quiver(q2u(1,:),q2u(2,:),gPhi3In(1,:)*scale,gPhi3In(2,:)*scale, 'Autoscale', 'off' ,'color','g')

sgtitle("gradiat of phi")

%N
figure(4), clf
N1 = phi1.*gPhi2-phi2.*gPhi1;
N2 = phi2.*gPhi3-phi3.*gPhi2;
N3 = phi3.*gPhi1-phi1.*gPhi3;

N1In = phi1In.*gPhi2In-phi2In.*gPhi1In;
N2In = phi2In.*gPhi3In-phi3In.*gPhi2In;
N3In = phi3In.*gPhi1In-phi1In.*gPhi3In;


subplot(3,1,1)
quiver(u,v,N1(1,:)*scale,N1(2,:)*scale,'Autoscale', 'off'), hold on
quiver(q2u(1,:),q2u(2,:),N1In(1,:)*scale,N1In(2,:)*scale,'Autoscale', 'off')
subplot(3,1,2)
quiver(u,v,N2(1,:)*scale,N2(2,:)*scale,'Autoscale', 'off'), hold on
quiver(q2u(1,:),q2u(2,:),N2In(1,:)*scale,N2In(2,:)*scale,'Autoscale', 'off')
subplot(3,1,3)
quiver(u,v,N3(1,:)*scale,N3(2,:)*scale,'Autoscale', 'off'), hold on
quiver(q2u(1,:),q2u(2,:),N3In(1,:)*scale,N3In(2,:)*scale,'Autoscale', 'off')

sgtitle("N functions")

%n x N (S)
figure(5), clf
S1 = [-(phi1.*gPhi2(2,:)-phi2.*gPhi1(2,:)); (phi1.*gPhi2(1,:)-phi2.*gPhi1(1,:))] ;
S2 = [-(phi2.*gPhi3(2,:)-phi3.*gPhi2(2,:)); (phi2.*gPhi3(1,:)-phi3.*gPhi2(1,:))] ; 
S3 = [-(phi3.*gPhi1(2,:)-phi1.*gPhi3(2,:)); (phi3.*gPhi1(1,:)-phi1.*gPhi3(1,:))] ;
col = [0.75,0.75,0.75];
subplot(3,1,1)
quiver(u,v,N1(1,:),N1(2,:),'color',col), hold on
quiver(u,v,S1(1,:),S1(2,:))
subplot(3,1,2)
quiver(u,v,N2(1,:),N2(2,:),'color',col), hold on
quiver(u,v,S2(1,:),S2(2,:))
subplot(3,1,3)
quiver(u,v,N3(1,:),N3(2,:),'color',col), hold on
quiver(u,v,S3(1,:),S3(2,:))


sgtitle("S functions")
% usnn{1} = [-up{2}*ug{1}(1)+up{1}*ug{2}(1); -up{2}*ug{1}(2)+up{1}*ug{2}(2); 0,0,0];
% usnn{2} = [-up{2}*ug{3}(1)+up{2}*ug{3}(1); -up{2}*ug{3}(2)+up{2}*ug{3}(2); 0,0,0];
% usnn{3} = [-up{3}*ug{1}(1)+up{1}*ug{3}(1); -up{3}*ug{1}(2)+up{1}*ug{3}(2); 0,0,0];
%(n x S)
figure(6), clf
nS1 = [-phi1.*gPhi2(1,:)+phi2.*gPhi1(1,:); -phi1.*gPhi2(2,:)+phi2.*gPhi1(2,:)] ;
nS2 = [-phi2.*gPhi3(1,:)+phi3.*gPhi2(1,:); -phi2.*gPhi3(2,:)+phi3.*gPhi2(2,:)] ; %fix in code
nS3 = [-phi3.*gPhi1(1,:)+phi1.*gPhi3(1,:); -phi3.*gPhi1(2,:)+phi1.*gPhi3(2,:)] ;
subplot(3,1,1)
quiver(u,v,S1(1,:),S1(2,:),'color',col), hold on
quiver(u,v,nS1(1,:),nS1(2,:))
subplot(3,1,2)
quiver(u,v,S2(1,:),S2(2,:),'color',col), hold on
quiver(u,v,nS2(1,:),nS2(2,:))
subplot(3,1,3)
quiver(u,v,S3(1,:),S3(2,:),'color',col), hold on
quiver(u,v,nS3(1,:),nS3(2,:))

sgtitle("n x S functions")


figure(7), clf, hold on
a = 10;
ranTri = rand(2,3)*a;
subplot(1,1,1)
plot([0,0,a,a,0],[0,a,a,0,0])
plot([ranTri(1,:),ranTri(1,1)],[ranTri(2,:),ranTri(2,1)])

scale = 0.4;
res = 0.5;
x = 0:res:a;
y = 0:res:a;
[x,y] = meshgrid(x,y);
q = zeros(size(y));
w = sin((pi*x)./a);
quiver(x,y,q*scale,w*scale,'Autoscale', 'off')

%q2u = [0,0;1,0;0,1]';

q2uTx = ranTri(1,1)+(ranTri(1,2)-ranTri(1,1))*q2u(2,:)+(ranTri(1,3)-ranTri(1,1))*q2u(1,:);
q2uTy = ranTri(2,1)+(ranTri(2,2)-ranTri(2,1))*q2u(2,:)+(ranTri(2,3)-ranTri(2,1))*q2u(1,:);
scatter(q2uTx,q2uTy);

qIn = zeros(size(q2uTx));
wIn = [1,1,1];%sin((pi*q2uTx)./a);

quiver(q2uTx,q2uTy,qIn*scale,wIn*scale,'g','Autoscale', 'off')

xlim([-1,11])
ylim([-1,11])

%%

% ipTmp = sum([qIn; wIn; qIn]);
% ipTmp * q2w' * det_jac;
%  
ug{1} = [-1 -1]';
ug{2} = [+1 0]';
ug{3} = [0 +1]';

xyz = [ranTri;0,0,0];
%Jacobian
jac = zeros(3);
for iIdx = 1:3
    jac = jac ...
        + [xyz(1,iIdx)*ug{iIdx}', 0; ...
           xyz(2,iIdx)*ug{iIdx}', 0; ...
           xyz(3,iIdx)*ug{iIdx}', 1]';
end

%hige res points in triangle

jacU = zeros(3);
jacU = [xyz(1,3)-xyz(1,3), xyz(1,2)-xyz(1,1); ...
        xyz(2,3)-xyz(2,1), xyz(2,2)-xyz(2,1)];

q2w = [1.666666666666667e-02; ...
       1.666666666666667e-02; ...
       1.666666666666667e-02]';
   
   
Tx = ranTri(1,1)+(ranTri(1,2)-ranTri(1,1))*u+(ranTri(1,3)-ranTri(1,1))*v;
Ty = ranTri(2,1)+(ranTri(2,2)-ranTri(2,1))*u+(ranTri(2,3)-ranTri(2,1))*v;

scatter(Tx,Ty,'x')

Area = 1/2*abs(xyz(1,1)*(xyz(2,2)-xyz(2,3))+ xyz(1,2)*(xyz(2,3)-xyz(2,1)) + xyz(1,3)*(xyz(2,1)-xyz(2,2)))
Area_estimate = [1,1,1]*q2w'*det(jacU)

% 
% num_of_points = 100;
% Area_pr_point = Area/100;
% 
% sum(points*Area_pr_point)
