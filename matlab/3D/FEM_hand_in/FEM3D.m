clear all
close all

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

direction = [1,-1];
a = 0.2;
b = 0.1;
c0 = 299792458; % speed of light in vacuum


% Materials not working
ma2er = {@(x,y,z) 1};
ma2si = {@(x,y,z) 0};

% Read mesh
file_list = ["waveguide_model"];
vers = 1;
load(sprintf('%s',file_list(vers)))
save_folder = 'test1';
%frequencies to simulate.
f_list = (0.1:0.4:1.5)*10^9; 

if ~exist(sprintf('/res/%s/%s',file_list(vers),save_folder), 'dir')
     disp('creates directory')
     mkdir(sprintf('res/%s/%s',file_list(vers),save_folder))
end

% Initialize the FEM
Fem_Init(no2xyz, ed2no_all, fa2no_all)
%ed2no_boundery = ed2no_pec;
%ed2no_pec = ed2no_bound;

% Find PEC edges in the database
edIdx_pec = ElementDatabase_Get('edges', ed2no_pec); % each edge where pec is present gets an id
noIdx_pec = unique(ed2no_pec(:))'; % each node where pec is present gets an id

% Find all edges in the database
edNum_all = ElementDatabase_Cardinal('edges'); % total number of edges
edIdx_all = 1:edNum_all; % each edge gets an id
noIdx_all = 1:size(no2xyz,2); % each node gets an id

%plot stucture with ports
 
% plot_edges(ed2no_port{1},no2xyz,13,'r');
% plot_edges(ed2no_port{2},no2xyz,13,'r');
% plot_edges(ed2no_port{3},no2xyz,13,'r');
% plot_edges(ed2no_pec,no2xyz,13,'g');

% Compute the interior edges
edIdx_int = setdiff(edIdx_all, edIdx_pec); % removes all edges that are pec from the index
noIdx_int = setdiff(noIdx_all, noIdx_pec); % removes all nodes that are pec from the index

%
S_par = zeros(length(f_list),length(port_fac2no_list),length(port_fac2no_list));
%parpool(2)
tic
%f_list = f_list(10:16);
listIdx = circshift(1:length(port_fac2no_list),-1);


for port_i = 1:1%length(port_fac2no_list)
    for fi = 1:length(f_list)
        tic
        f = f_list(fi);
        fprintf('frequency: %.4f GHz\n',f*10^(-9))

        w = f*2*pi;
        k0 = w*(1/c0);
        k_z10 = sqrt((k0^2)-(pi/a)^2);
        Z = (w*(pi*4*10^(-7))/k_z10);
        watt = 0.5;
        E0 = sqrt(watt*Z);
        gamma = 1j*k_z10;
        

        [KeMtx, BeMtx, bMtx] = ...
            Fem_Assemble(no2xyz, el2no, el2ma, ma2er, ma2si, port_fac2no_list, k0, gamma, k_z10, a, E0, direction);

        KeMtx = KeMtx(edIdx_int,edIdx_int);

        bMtx = bMtx(edIdx_int);
        BeMtx = BeMtx(edIdx_int,edIdx_int);
        KMtx = KeMtx+BeMtx;

        E = KMtx\bMtx;
        
        eFld_all = zeros(edNum_all,1); % prealocates memmory and includes edges which are PEC. PEC are zero
        eFld_all(edIdx_int) = E;

        S_par(fi,:,port_i) = S_parameters(eFld_all,port_fac2no_list,no2xyz,a,b,k_z10,E0);

        fprintf("S11: %f + %fi \tS12: %f +%fi \n",real(S_par(fi,1,port_i)),imag(S_par(fi,1,port_i)),real(S_par(fi,2,port_i)),imag(S_par(fi,2,port_i)))
        filename = sprintf('res/%s/%s/E_filds_%d_f_%.0f',file_list(vers),save_folder,port_i,f*10^-6);
        %save(filename,'eFld_all','ed2no_boundery','no2xyz','el2no');
        save(filename,'eFld_all','ed2no_pec','ed2no_port','no2xyz','el2no');
        toc
    end
    direction = direction(listIdx);
    port_fac2no_list = port_fac2no_list(listIdx);
end

save(sprintf('res/%s/%s/Sparamters',file_list(vers),save_folder),'S_par','f_list')

figure(4), clf
subplot(2,1,1)
plot(f_list*10^(-9),10*log10(abs(S_par(:,1,1))),'DisplayName','S11')
hold on
plot(f_list*10^(-9),10*log10(abs(S_par(:,2,1))),'DisplayName','S21')
xlabel('Freqency [GHz]')
ylabel('dB')
lg = legend();
lg.Location = 'northwest';
ylim([-60,1])

subplot(2,1,2)
plot(f_list*10^(-9),(angle(S_par(:,1,1)))*180/pi,'DisplayName','S11')
hold on
plot(f_list*10^(-9),(angle(S_par(:,2,1)))*180/pi,'DisplayName','S21')
xlabel('Freqency [GHz]')
ylabel('Phase [$\degree$]')
set(gcf,'Position',[100 100 600 600]);


saveas(gcf,'res/Sparameters_our_FEM','svg')




