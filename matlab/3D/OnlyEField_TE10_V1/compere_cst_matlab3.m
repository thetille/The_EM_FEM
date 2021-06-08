clear
close all

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

file_list = ["cylinder_waveguide2", "waveguide_model3 - simple"...
,"waveguide_model3_flat7","mesh_cylinder_R0"...
,"waveguide_model3_flathigh","waveguide_model3_wired"...
,"waveguide_model3_highHigh","waveguide_model3_flat_long"];
vers = 3;
save_folder = 'test1';
port_i = 1;

matlab_mesh = load(file_list(vers));
Fem_Init(matlab_mesh.no2xyz, matlab_mesh.ed2no_all, matlab_mesh.fa2no_all);
pMtx_ed2no = ProjSol2Nodes_Assemble(matlab_mesh.no2xyz, matlab_mesh.el2no);

max1 = 0; max2 = 0; max3 = 0; max4 = 0;
min1 = 1000; min2 = 1000; min3 = 1000; min4 = 1000;

f_list = (0.25:0.30:1.45)*10^9;
font_size = 11.5;
first = true;
i = 0;
normelize = load('normelization.mat');

for f = f_list%0.75*10^9
    if f < 0.75*10^9
        figure_offset = 0;
        cut_off = 'cut';
    else
       figure_offset = 2; 
       cut_off = '';
       if first == true
           %max1 = 0; max2 = 0; max3 = 0; max4 = 0;
           %min1 = 1000; min2 = 1000; min3 = 1000; min4 = 1000;
           first = false;
       end
    end
    
    i = i +1;
    %sgtitle(sprintf('f: %f GHz',f*10^-9))
    
    %%%%% Analytical solution %%%%%%
    res = 0.01;
    m = 1;
    n = 0;
    omega = f*2*pi;
    mu = 4*pi*10^(-7);
    eps = 8.8541878128*10^(-12);
    b = 0.1;
    a = 0.2;
    kc = sqrt(((m*pi/a)^2)+((n*pi/b)^2));
    k = omega*sqrt(mu*eps);
    be = sqrt((k^2) - (kc^2));
    
%     eta = sqrt(mu/eps);
%     Z = ((k*eta)/be);
%     E0 = sqrt(0.5*Z);
    
    A = (sqrt(2)*pi)/((a^(2/3))*sqrt(b)*sqrt(mu)*sqrt(omega)*sqrt(be));
    
    %(omega*mu*(a^3)*(A^2)*b)/(
    
    x = 0.05; %the anlytical has the (0,0,0) point at the center inted of the middle
    y = 0.025; %the anlytical has the (0,0,0) point at the center inted of the middle
    z = 0:res:0.4;    
    Ex_ana = 1i*(( 1j*omega*mu*n*pi)/((kc^2)*b)) *A* cos((m*pi*x)/a) * sin((n*pi*y)/b) * exp(1j*be.*z);
    Ey_ana = 1i*((-1j*omega*mu*m*pi)/((kc^2)*a)) *A* sin((m*pi*x)/a) * cos((n*pi*y)/b) * exp(1j*be.*z);
    Ez_ana = zeros(size(Ey_ana));
    
    
    %%%%%% loade CST FEM data %%%%%
    if mod(f*10^-9,1) < 0.0001
        filename = sprintf('test/cst_res/e-field (f=%.0f) [1].txt',f*10^-9);
    elseif mod(f*10^-9,0.1) > 0.01
        filename = sprintf('test/cst_res/e-field (f=%.2f) [1].txt',f*10^-9);
    else
        filename = sprintf('test/cst_res/e-field (f=%.1f) [1].txt',f*10^-9);
    end
    
    delimiterIn = ' ';
    headerlinesIn = 2;
    CST_data = importdata(filename,delimiterIn,headerlinesIn);
    CST_data = CST_data.data; %x,y,z, ExRe,EyRe,EzRe, ExIm,EyIm,EzIm 
    

    CST_Exabs = abs(CST_data(:,4)+CST_data(:,4+3)*1i);
    CST_Exphase = angle(CST_data(:,4)+CST_data(:,4+3)*1i);
    
    CST_Eyabs = abs(CST_data(:,5)+CST_data(:,5+3)*1i);
    CST_Eyphase = angle(CST_data(:,5)+CST_data(:,5+3)*1i);
    
    CST_Ezabs = abs(CST_data(:,6)+CST_data(:,6+3)*1i);
    CST_Ezphase = angle(CST_data(:,6)+CST_data(:,6+3)*1i);

    CSTz = CST_data(:,3);
    
    %%%%%% loade CST FIT data %%%%%
    mode = 1;
    if mod(f*10^-9,1) < 0.0001
        filename = sprintf('test/cst_res_fit/e-field (f=%.0f) [1].txt',f*10^-9);
    elseif mod(f*10^-9,0.1) > 0.01
       filename = sprintf('test/cst_res_fit/e-field (f=%.2f) [1].txt',f*10^-9);
    else
       filename = sprintf('test/cst_res_fit/e-field (f=%.1f) [1].txt',f*10^-9);
    end
    delimiterIn = ' ';
    headerlinesIn = 2;
    CST_data_fit = importdata(filename,delimiterIn,headerlinesIn);
    CST_data_fit = CST_data_fit.data; %x,y,z, ExRe,EyRe,EzRe, ExIm,EyIm,EzIm 
    

    CST_Exabs_fit = abs(CST_data_fit(:,4)+CST_data_fit(:,4+3)*1i);
    CST_Exphase_fit = angle(CST_data_fit(:,4)+CST_data_fit(:,4+3)*1i);
    
    CST_Eyabs_fit = abs(CST_data_fit(:,5)+CST_data_fit(:,5+3)*1i);
    CST_Eyphase_fit = angle(CST_data_fit(:,5)+CST_data_fit(:,5+3)*1i);
    
    CST_Ezabs_fit = abs(CST_data_fit(:,6)+CST_data_fit(:,6+3)*1i);
    CST_Ezphase_fit = angle(CST_data_fit(:,6)+CST_data_fit(:,6+3)*1i);

    CSTz_fit = CST_data_fit(:,3);
    
    %%%%% load matlab data %%%%%
    filename = sprintf('res/%s/%s/E_filds_%d_f_%.0f.mat',file_list(vers),save_folder,port_i,f*10^-6);
    matlab_data = load(filename);
    
    %%% interpolate Matlba data %%%%%
    
    %res2 = -res;
    Zq = min(matlab_data.no2xyz(3,:)):res:max(matlab_data.no2xyz(3,:));
    
    Yq = ones(size(Zq))*0.025;
    Xq = ones(size(Zq))*-0.05;
    
    exFld_all = pMtx_ed2no.xc*matlab_data.eFld_all;
    eyFld_all = pMtx_ed2no.yc*matlab_data.eFld_all;
    ezFld_all = pMtx_ed2no.zc*matlab_data.eFld_all;
    
    matlab_ExAbs = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),abs(exFld_all*(-1i)),Xq,Yq,Zq);
    matlab_ExPhase = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),angle(exFld_all*(-1i)),Xq,Yq,Zq);
    
    matlab_EyAbs = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),abs(eyFld_all*(-1i)),Xq,Yq,Zq);
    matlab_EyPhase = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),angle(eyFld_all*(-1i)),Xq,Yq,Zq);
    
    matlab_EzAbs = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),abs(ezFld_all*(-1i)),Xq,Yq,Zq);
    matlab_EzPhase = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),angle(ezFld_all*(-1i)),Xq,Yq,Zq);
    
    ax = subplot(3,2,i); hold on
    title(sprintf('Frequency: %.2f',f*10^-9))
    %plot(z,abs(Ey_ana),'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    plot(CSTz,CST_Eyabs*normelize.Eyabs_ana_CSTFEM(i),'DisplayName','CST FEM')
    plot(CSTz_fit,CST_Eyabs_fit*normelize.Eyabs_ana_CSTFIT(i),'DisplayName','CST FIT')
    plot(Zq,matlab_EyAbs*normelize.Eyabs_ana_mat(i),'DisplayName','Our FEM')
    xlabel('Z cordinate [m]')
    ylabel('E-Field [V/m]')
   
    fprintf('var f:%.2f,CST FEM: %f\n',f*10^-9,var(CST_Eyabs*normelize.Eyabs_ana_CSTFEM(i)))
    fprintf('var f:%.2f,CST FIT: %f\n',f*10^-9,var(CST_Eyabs_fit*normelize.Eyabs_ana_CSTFIT(i)))
    fprintf('var f:%.2f,Our FEM: %f\n',f*10^-9,var(matlab_EyAbs*normelize.Eyabs_ana_mat(i)))
    
end


lh =legend(ax,'Location','SouthOutside');%'Orientation','Horizontal');
lh.FontSize = 11.5;