clear all
close all

global normals;
normals = 0; %enable plotting of port normals
% Direct or sparse eigenvalue solver for small and large
% problems, respectively [solver = 'direct' or 'sparse']

a = 0.2;
b = 0.1;
c0 = 299792458; % speed of light in vacuum


% Materials
ma2er = {@(x,y,z) 1};% + 4*exp(-((x-0.125).^2+y.^2+(z-0.3).^2)/(0.1^2))};
ma2si = {@(x,y,z) 0};%0.1*exp(-((x+0.125).^2+y.^2+(z-0.3).^2)/(0.1^2))};

% Read mesh
file_list = ["cylinder_waveguide2", "waveguide_model3 - simple"...
            ,"waveguide_model3_flat","mesh_cylinder_R0"...
            ,"waveguide_model3_flathigh","waveguide_model3_wired"...
            ,"waveguide_model3_highHigh","waveguide_model3_flat_long"];
vers = 3;
load(file_list(vers))
save_folder = 'test5';

if ~exist(sprintf('/res/%s/%s',file_list(vers),save_folder), 'dir')
     disp('creates directory')
    mkdir(sprintf('res/%s/%s',file_list(vers),save_folder))
end

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
edIdx_all = 1:edNum_all; % each edge gets an id
noIdx_all = 1:size(no2xyz,2); % each node gets an id

%plot stucture with ports
 
% plot_edges(ed2no_port1,no2xyz,13,'r');
% plot_edges(ed2no_port2,no2xyz,13,'g');
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

f_list = (0.6:0.025:1.2)*10^9;%[0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]*10^9;
%f_list = 0.6*10^9;
%f_list = 1.2*10^9;
%
S_par = zeros(length(f_list),2);
%parpool(2)
tic
%f_list = f_list(10:16);
for fi = 1:length(f_list)
    f = f_list(fi);
    fprintf('frequency: %.4f GHz\n',f*10^(-9))
    
    w = f*2*pi;
    k0 = w*(1/c0);
    k_z10 = sqrt((k0^2)-(pi/a)^2);
    Z = (w*(pi*4*10^(-7))/k_z10)
    E0 = sqrt(0.5*Z)
    gamma = 1j*k_z10;
    
%     if (pi/a)^2 <= k0^2
%         %k_z10 = sqrt(((pi/a).^2)-k0.^2);
%         gamma = sqrt(((pi/a).^2)-k0.^2);
%         disp('gamma2')
%     else
%         %k_z10 = sqrt(k0^2-(pi/a)^2);
%         gamma = 1j*k_z10;
%         disp('gamma1')
%     end

    [KeMtx, BeMtx, bMtx] = ...
        Fem_Assemble(no2xyz, el2no, el2ma, ma2er, ma2si, fac2no_port1, fac2no_port2, k0, gamma, k_z10, a, E0);

    KeMtx = KeMtx(edIdx_int,edIdx_int);
    
    bMtx = bMtx(edIdx_int);
    BeMtx = BeMtx(edIdx_int,edIdx_int);
    %KeMtx(edIdx_port1_int,edIdx_port1_int) = 1;
    %KeMtx(edIdx_port2_int,edIdx_port2_int) = 0;
    KMtx = KeMtx+BeMtx;
    

    tic
    E = KMtx\bMtx;
    toc

    eFld_all = zeros(edNum_all,1); % prealocates memmory and includes edges which are PEC. PEC are zero
    eFld_all(edIdx_int) = E;%diag(KMtx);%diag(BeMtx);
    %eFld_all(edIdx_int) = bMtx;


    S_par(fi,:) = S_parameters(eFld_all,fac2no_port1,fac2no_port2,no2xyz,a,b,k_z10,E0);
    
    %fprintf("S11: %f \tS12: %f \n",abs(S_par(fi,1)),abs(S_par(fi,2)))
    fprintf("S11: %f + %fi \tS12: %f \n",real(S_par(fi,1)),imag(S_par(fi,1)),abs(S_par(fi,2)))
    filename = sprintf('res/%s/%s/E_filds_f_%.0f',file_list(vers),save_folder,f*10^-6);
    save(filename,'eFld_all','ed2no_boundery','no2xyz','el2no');
    %plot_fileds(eFld_all,ed2no_boundery,no2xyz)
end
toc
figure(1), clf
plot(f_list*10^(-9),abs(S_par(:,1)),'DisplayName','S11')
hold on
plot(f_list*10^(-9),abs(S_par(:,2)),'DisplayName','S12')
legend()
%ylim([0,2])
save(sprintf('res/%s/%s/Sparamters',file_list(vers),save_folder),'S_par','f_list')

