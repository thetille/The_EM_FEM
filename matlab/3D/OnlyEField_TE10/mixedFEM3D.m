clear all
figure(10), clf
global normals;
normals = 0; %enable plotting of port normals
% Direct or sparse eigenvalue solver for small and large
% problems, respectively [solver = 'direct' or 'sparse']

a = 0.2;
c0 = 299792458;
f = 1*10^9; % 1 ghz
w = f*2*pi;
k0 = w*(1/c0);
k_z10 = sqrt(k0^2-(pi/a)^2);
gamma = 1j*k_z10;


% Materials
ma2er = {@(x,y,z) 1};% + 4*exp(-((x-0.125).^2+y.^2+(z-0.3).^2)/(0.1^2))};
ma2si = {@(x,y,z) 0};%0.1*exp(-((x+0.125).^2+y.^2+(z-0.3).^2)/(0.1^2))};

% Constants
% c0 = 299792458;     % speed of light in vacuum
% m0 = 4*pi*1e-7;     % permeability in vacuum
% e0 = 1/(m0*c0^2);   % permittivity in vacuum
% z0 = sqrt(m0/e0);   % wave 	 in vacuum

% Read mesh
file_list = ["cylinder_waveguide2", "waveguide_model3 - simple","waveguide_model3","mesh_cylinder_R0","waveguide_model3_highres"];
load(file_list(3))

% ed2no_pec = [ed2no_port1, ed2no_port2, ed2no_bound];

% Initialize the FEM
Fem_Init(no2xyz, ed2no_all, fa2no_all)
ed2no_boundery = ed2no_pec;
ed2no_pec = ed2no_bound;
% Find PEC edges in the database
edIdx_pec = ElementDatabase_Get('edges', ed2no_pec); % each edge where pec is present gets an id
edIdx_port1 = ElementDatabase_Get('edges', ed2no_port1);
edIdx_port2 = ElementDatabase_Get('edges', ed2no_port2);
noIdx_pec = unique(ed2no_pec(:))'; % each node where pec is present gets an id

% Find all edges in the database
edNum_all = ElementDatabase_Cardinal('edges'); % total number of edges
%faNum_all = ElementDatabase_Cardinal('faces'); % total number of faces
%only for H field calculations
edIdx_all = 1:edNum_all; % each edge gets an id
noIdx_all = 1:size(no2xyz,2); % each node gets an id

%plot stucture with ports
% plot_edges(ed2no_port1,no2xyz,13,'r');
% plot_edges(ed2no_port2,no2xyz,13,'r');
% plot_edges(ed2no_pec,no2xyz,13,'g');

% Compute the interior edges
edIdx_int = setdiff(edIdx_all, edIdx_pec); % removes all edges that are pec from the index
noIdx_int = setdiff(noIdx_all, noIdx_pec); % removes all nodes that are pec from the index

%varibles for port edges
edIdx_port1_int = setdiff(edIdx_port1,edIdx_pec);
edIdx_port2_int = setdiff(edIdx_port2,edIdx_pec);
% Assemble global matrices

%plot to show results, needs to be here in order for normals debug code to
%work
figure(3), clf;
subplot(1,2,1), hold on

[KeMtx, BeMtx, bMtx] = ...
    Fem_Assemble(no2xyz, el2no, el2ma, ma2er, ma2si, fac2no_port1, fac2no_port2, k0, gamma, k_z10, a);

% no2xyz = coordinates to all points
% el2no = all points in a tetra
% el2ma = material indices of the tetrahedrons
% ma2er = stores the permittivity associated with the different material indices
% ma2si = storesthe conductivity associated with the different material indices
% fac2no_port1 = The triangles in port1
% fac2no_port2 = The triangles in port2

% sparce to full
%KeMtx = full(KeMtx);
%BeMtx = full(BeMtx);

%set pec bounderys (edge elements) to 0
% KeMtx(edIdx_pec,edIdx_pec) = 0;
% KeMtx(edIdx_port1,edIdx_port1) = 1;
% KeMtx(edIdx_port2,edIdx_port2) = 0;

%KeMtx(edIdx_port1_int,edIdx_port1_int) = 0;
%KeMtx(edIdx_port2_int,edIdx_port2_int) = 0;
KeMtx = KeMtx(edIdx_int,edIdx_int);


%bMtx(intersect(edIdx_port1,edIdx_pec)) = 0;
bMtx = bMtx(edIdx_int);
% bMtx = diag(bMtx);
%BeMtx_inter = intersect([edIdx_port1,edIdx_port2],edIdx_pec);
%BeMtx(BeMtx_inter,BeMtx_inter) = 0;
%edIdx_int = setdiff([edIdx_port1,edIdx_port2],edIdx_pec);
BeMtx = BeMtx(edIdx_int,edIdx_int);
KMtx = KeMtx+BeMtx;

pMtx_ed2no = ProjSol2Nodes_Assemble(no2xyz, el2no);

%%
%A = null(KMtx);
%E = A(:,1);
E = KMtx\bMtx;

eFld_all = zeros(edNum_all,1); % prealocates memmory and includes edges which are PEC. PEC are zero
eFld_all(edIdx_int) = E;
%eFld_all(edIdx_int) = bMtx;

%E(edIdx_pec) = 0;
exFld_all = pMtx_ed2no.xc*eFld_all;
eyFld_all = pMtx_ed2no.yc*eFld_all;
ezFld_all = pMtx_ed2no.zc*eFld_all;
% 
% exFld_all = pMtx_ed2no.xc*bMtx';
% eyFld_all = pMtx_ed2no.yc*bMtx';
% ezFld_all = pMtx_ed2no.zc*bMtx';
%dVal = max(no2xyz,[],'all')*1.2;

%get fields at port 1
noIdx_port1 = unique(fac2no_port1(:));
exFld_port1 = exFld_all(noIdx_port1);
eyFld_port1 = eyFld_all(noIdx_port1);
ezFld_port1 = ezFld_all(noIdx_port1);

%get fields at port 2
noIdx_port2 = unique(fac2no_port2(:));
exFld_port2 = exFld_all(noIdx_port2);
eyFld_port2 = eyFld_all(noIdx_port2);
ezFld_port2 = ezFld_all(noIdx_port2);

S_parameters(E,ed2no_port1,ed2no_port2,fac2no_port1,fac2no_port2,no2xyz,a)

figure(3)
subplot(1,2,1)

for edIdx = 1:size(ed2no_boundery,2)
    noTmp = ed2no_boundery(:,edIdx);
    xyzTmp = no2xyz(:,noTmp);
    plot3(xyzTmp(1,:), xyzTmp(2,:), xyzTmp(3,:), ...
        'Color', 0.5*[1 1 1])
end
exViz = real(exFld_all(:).');
eyViz = real(eyFld_all(:).');
ezViz = real(ezFld_all(:).');
quiverC3D(no2xyz(1,:)', no2xyz(2,:)', no2xyz(3,:)', ...
    exViz', eyViz', ezViz',1)
set(gca,'Color',[0,27/100,55/100])
axis equal
%axis off
view(40,-20)

subplot(1,2,2), hold on

for edIdx = 1:size(ed2no_boundery,2)
    noTmp = ed2no_boundery(:,edIdx);
    xyzTmp = no2xyz(:,noTmp);
    plot3(xyzTmp(1,:), xyzTmp(2,:), xyzTmp(3,:), ...
        'Color', 0.5*[1 1 1])
end
exViz = imag(exFld_all(:).');
eyViz = imag(eyFld_all(:).');
ezViz = imag(ezFld_all(:).');
quiverC3D(no2xyz(1,:)', no2xyz(2,:)', no2xyz(3,:)', ...
    exViz', eyViz', ezViz',1)
set(gca,'Color',[0,27/100,55/100])

axis equal
%axis off
view(40,-20)

X =  no2xyz(1,:);
Y =  no2xyz(2,:);
Z =  no2xyz(3,:);

V = abs(exViz)+abs(eyViz)+abs(ezViz);

res = [0.005,0.005,0.005];
[Xq,Yq,Zq] = meshgrid(-0.1:res(1):0.1, -.1:res(2):.1 ,0:res(3):0.4);
Vq = griddata(X,Y,Z,V,Xq,Yq,Zq);

figure(5), clf;
slice(Xq,Yq,Zq,Vq,[-0.1 0.1],[-0.1 0.1],[0 0.4]);
shading flat
pbaspect([1 1 1])
view(40,-14+(rand(1)*0.1))
