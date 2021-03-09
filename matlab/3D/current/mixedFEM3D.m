clear all

% Direct or sparse eigenvalue solver for small and large
% problems, respectively [solver = 'direct' or 'sparse']
solver = 'direct';

% Materials
ma2er = {@(x,y,z) 1 + 4*exp(-((x-0.125).^2+y.^2+(z-0.3).^2)/(0.1^2))};
ma2si = {@(x,y,z) 0.1*exp(-((x+0.125).^2+y.^2+(z-0.3).^2)/(0.1^2))};

% Constants
c0 = 299792458;     % speed of light in vacuum
m0 = 4*pi*1e-7;     % permeability in vacuum
e0 = 1/(m0*c0^2);   % permittivity in vacuum
z0 = sqrt(m0/e0);   % wave impedance in vacuum

% Read mesh
load mesh_cylinder_R0

% Initialize the FEM
Fem_Init(no2xyz, ed2no_all, fa2no_all)

% Find PEC edges in the database
edIdx_pec = ElementDatabase_Get('edges', ed2no_pec); % each edge where pec is present gets an id
noIdx_pec = unique(ed2no_pec(:))'; % each node where pec is present gets an id

% Find all edges in the database
edNum_all = ElementDatabase_Cardinal('edges'); % total number of edges
faNum_all = ElementDatabase_Cardinal('faces'); % total number of faces
edIdx_all = 1:edNum_all; % each edge gets an id
noIdx_all = 1:size(no2xyz,2); % each node gets an id

% Compute the interior edges
edIdx_int = setdiff(edIdx_all, edIdx_pec); % removes all edges that are pec from the index
noIdx_int = setdiff(noIdx_all, noIdx_pec); % removes all nodes that are pec from the index

% Assemble global matrices
[eMtx, sMtx, cMtx, uMtx] = ...
    Fem_Assemble(no2xyz, el2no, el2ma, ma2er, ma2si);

% no2xyz = coordinates to all points
% el2no = all points in a tetra
% el2ma = material indices of the tetrahedrons
% ma2er = stores the permittivity associated with the different material indices
% ma2si = storesthe conductivity associated with the different material indices

eMtx_int = eMtx(edIdx_int,edIdx_int);    % mass matrix with permittivity
sMtx_int = sMtx(edIdx_int,edIdx_int);    % mass matrix with conductivity
cMtx_int = cMtx(:,edIdx_int);            % curl matrix
uMtx_int = uMtx;                         % mass matrix with unity coefficient

edNum_dof = length(edIdx_int);
faNum_dof = faNum_all;

solIdx_bFld = 1:faNum_dof;
solIdx_eFld = faNum_dof + (1:edNum_dof);

aMtx = ...
    [sparse(faNum_dof,faNum_dof) -cMtx_int; ...
    cMtx_int.' -z0*sMtx_int];
bMtx = ...
    [uMtx_int sparse(faNum_dof,edNum_dof); ...
    sparse(edNum_dof,faNum_dof) eMtx_int];

% Solve the eigenvalue problem
if strcmp(solver, 'direct')
    aMtx = full(aMtx);
    bMtx = full(bMtx);
    [eigVtr_int, eigVal] = eig(aMtx, bMtx);
elseif strcmp(solver, 'sparse')
    [eigVtr_int, eigVal] = eigs(aMtx, 0.5*(bMtx+bMtx'), 30, 1i*(19.5+1i));
else
    error('unknown eigenvalue solver')
end

eigVal = diag(eigVal); % j*w/c0
[eigTmp, eigIdx_sort] = sort(real(-1i*eigVal));
eigVal     = -1i*eigVal(eigIdx_sort); % w/c0
eigVtr_int = eigVtr_int(:,eigIdx_sort);
fr = c0*eigVal/(2*pi);


%% Dsiplays all data

% Select the modes that will be visualized
mVtr = 1151 + (1:20);

% Visualize the eigenfrequencies
figure(1), clf
plot(real(fr(mVtr))/1e9, imag(fr(mVtr))/1e9, 'ks')
xlabel('Real part of eigenfrequency [GHz]')
ylabel('Imaginary part of eigenfrequency [GHz]')
grid on

% Visualize the eigenmodes
bFld_all = eigVtr_int(solIdx_bFld,:);

eFld_all = zeros(edNum_all,size(eigVtr_int,2));
eFld_all(edIdx_int,:) = eigVtr_int(solIdx_eFld,:);

[pMtx_ed2no, pMtx_fa2no] = ProjSol2Nodes_Assemble(no2xyz, el2no);

bxFld_all = (pMtx_fa2no.xc*bFld_all) / c0;
byFld_all = (pMtx_fa2no.yc*bFld_all) / c0;
bzFld_all = (pMtx_fa2no.zc*bFld_all) / c0;

exFld_all = pMtx_ed2no.xc*eFld_all;
eyFld_all = pMtx_ed2no.yc*eFld_all;
ezFld_all = pMtx_ed2no.zc*eFld_all;



for mIdx = mVtr
    figure(2), clf
    dVal = 0.4;
    for dIdx = 0:1
        for edIdx = 1:size(ed2no_pec,2)
            noTmp = ed2no_pec(:,edIdx);
            xyzTmp = no2xyz(:,noTmp);
            plot3(dIdx*dVal + xyzTmp(1,:), xyzTmp(2,:), xyzTmp(3,:), ...
                'Color', 0.5*[1 1 1]), hold on
        end
        if dIdx == 0
            exViz = real(exFld_all(:,mIdx).');
            eyViz = real(eyFld_all(:,mIdx).');
            ezViz = real(ezFld_all(:,mIdx).');
            quiver3(dIdx*dVal + no2xyz(1,:), no2xyz(2,:), no2xyz(3,:), ...
                exViz, eyViz, ezViz, 2, 'k')
        elseif dIdx == 1
            bxViz = imag(bxFld_all(:,mIdx).');
            byViz = imag(byFld_all(:,mIdx).');
            bzViz = imag(bzFld_all(:,mIdx).');
            quiver3(dIdx*dVal + no2xyz(1,:), no2xyz(2,:), no2xyz(3,:), ...
                bxViz, byViz, bzViz, 2, 'k')
        end
    end
    axis equal
    axis off
    view(-24,14)
end