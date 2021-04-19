clear all

% Direct or sparse eigenvalue solver for small and large
% problems, respectively [solver = 'direct' or 'sparse']
solver = ['direct';'sparse'];
solver = solver(2,:);

% Materials
ma2er = {@(x,y,z) 1};% + 4*exp(-((x-0.125).^2+y.^2+(z-0.3).^2)/(0.1^2))};
ma2si = {@(x,y,z) 0};%0.1*exp(-((x+0.125).^2+y.^2+(z-0.3).^2)/(0.1^2))};

% Constants
c0 = 299792458;     % speed of light in vacuum
m0 = 4*pi*1e-7;     % permeability in vacuum
e0 = 1/(m0*c0^2);   % permittivity in vacuum
z0 = sqrt(m0/e0);   % wave 	 in vacuum

% Read mesh
file_list = ["cylinder_waveguide2", "waveguide_model3","mesh_cylinder_R0"];
load(file_list(2))

% ed2no_pec = [ed2no_port1, ed2no_port2, ed2no_bound];

% Initialize the FEM
Fem_Init(no2xyz, ed2no_all, fa2no_all)

% Find PEC edges in the database
edIdx_pec = ElementDatabase_Get('edges', ed2no_pec); % each edge where pec is present gets an id
edIdx_port1 = ElementDatabase_Get('edges', ed2no_port1);
edIdx_port2 = ElementDatabase_Get('edges', ed2no_port2);
noIdx_pec = unique(ed2no_pec(:))'; % each node where pec is present gets an id

% Find all edges in the database
edNum_all = ElementDatabase_Cardinal('edges'); % total number of edges
faNum_all = ElementDatabase_Cardinal('faces'); % total number of faces
edIdx_all = 1:edNum_all; % each edge gets an id
noIdx_all = 1:size(no2xyz,2); % each node gets an id

% Compute the interior edges
edIdx_int = setdiff(edIdx_all, [edIdx_pec,edIdx_port1,edIdx_port2]); % removes all edges that are pec from the index
noIdx_int = setdiff(noIdx_all, noIdx_pec); % removes all nodes that are pec from the index
tic
% Assemble global matrices
[KeMtx, BeMtx, bMtx] = ...
    Fem_Assemble(no2xyz, el2no, el2ma, ma2er, ma2si, fac2no_port1, fac2no_port2);

% no2xyz = coordinates to all points
% el2no = all points in a tetra
% el2ma = material indices of the tetrahedrons
% ma2er = stores the permittivity associated with the different material indices
% ma2si = storesthe conductivity associated with the different material indices

%set all bounderys to 0
KeMtx(edIdx_pec) = 0;
KeMtx(edIdx_port1) = 0;
KeMtx(edIdx_port2) = 0;

KMtx = KeMtx + BeMtx;
%E = bMtx\KMtx;



%% Dsiplays all data


% % Select the modes that will be visualized
% if solver == 'sparse'
%     mVtr = 1:30;
% else
%   mVtr = length(noIdx_int)+1:length(noIdx_int)+1+30;%length(fr)%length(no2xyz):length(fr);
% end

% 
% % Visualize the eigenfrequencies
% figure(1), clf
% plot(real(fr(mVtr))/1e9, imag(fr(mVtr))/1e9, 'ks')
% xlabel('Real part of eigenfrequency [GHz]')
% ylabel('Imaginary part of eigenfrequency [GHz]')
% grid on
% 
% figure(2), clf
% %subplot(2,1,1)
% plot(real(fr(mVtr))/1e9,'X','DisplayName','Implementation')
% %subplot(2,1,2)
% hold on
% load("cst_data")
% plot(cst_data,'X','DisplayName','CST')
% xlabel('Modes []')
% ylabel('Frequency [GHz]')
% legend
% grid on

% Visualize the eigenmodes
%solIdx_bFld = 1:faNum_dof; % The solution for the magnetic field  
%solIdx_eFld = faNum_dof + (1:edNum_dof); % the solution for the electric field

%bFld_all = eigVtr_int(solIdx_bFld,:); % interiour magnetic field

% eFld_all = zeros(edNum_all,size(E,2)); % prealocates memmory and includes edges which are PEC. PEC are zero
% eFld_all(edIdx_int,:) = E(:,:); % sets all interieor edges to eigenvector value

pMtx_ed2no = ProjSol2Nodes_Assemble(no2xyz, el2no);

%bxFld_all = (pMtx_fa2no.xc*bFld_all) / c0;
%byFld_all = (pMtx_fa2no.yc*bFld_all) / c0;
%bzFld_all = (pMtx_fa2no.zc*bFld_all) / c0;

exFld_all = pMtx_ed2no.xc*bMtx;
eyFld_all = pMtx_ed2no.yc*bMtx;
ezFld_all = pMtx_ed2no.zc*bMtx;

%%
%dVal = max(no2xyz,[],'all')*1.2;

figure(3), clf;

%for dIdx = 0:1
for edIdx = 1:size(ed2no_pec,2)
    noTmp = ed2no_pec(:,edIdx);
    xyzTmp = no2xyz(:,noTmp);
    plot3(xyzTmp(1,:), xyzTmp(2,:), xyzTmp(3,:), ...
        'Color', 0.5*[1 1 1]), hold on
end

    exViz = real(exFld_all(:).');
    eyViz = real(eyFld_all(:).');
    ezViz = real(ezFld_all(:).');
    quiver3(no2xyz(1,:), no2xyz(2,:), no2xyz(3,:), ...
        exViz, eyViz, ezViz, 2, 'k')

axis equal
axis off
view(-24,14)
%title( sprintf('eigenfrequency [GHz]: %2.2f + %2.2fi',real(fr(mIdx))/1e9,imag(fr(mIdx))/1e9 ))

exViz = real(exFld_all(:).');
eyViz = real(eyFld_all(:).');
ezViz = real(ezFld_all(:).');

X =  no2xyz(1,:);
Y =  no2xyz(2,:);
Z =  no2xyz(3,:);

% %phase_list = [1:-0.1:-1];
% %i = 1;
% %for phase = [phase_list flip(phase_list(1:end-1)) ]
%     figure(4), clf;
%     %V = (phase*exViz).^2+(phase*eyViz).^2+(phase*ezViz).^2;
%     V = (exViz).^2+(eyViz).^2+(ezViz).^2;
% 
%     res = [0.005,0.005,0.005];
%     [Xq,Yq,Zq] = meshgrid(-0.1:res(1):0.1, -.1:res(2):.1 ,0:res(3):0.4);
%     Vq = griddata(X,Y,Z,V,Xq,Yq,Zq);
% 
%     slice(Xq,Yq,Zq,Vq,0,0,[0.1 0.3]);
%     shading flat
%     title(sprintf('eigenfrequency [GHz]: %2.2f + %2.2fi',real(fr(mIdx))/1e9,imag(fr(mIdx))/1e9 ))
%     %F(i) = getframe();
%     %i = i+1;
% %end
% %figure(5)
% %movie(F,10)


