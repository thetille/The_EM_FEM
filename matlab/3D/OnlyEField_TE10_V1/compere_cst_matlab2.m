clear
close all

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

file_list = ["cylinder_waveguide2", "waveguide_model3 - simple"...
,"waveguide_model3_flat5","mesh_cylinder_R0"...
,"waveguide_model3_flathigh","waveguide_model3_wired"...
,"waveguide_model3_highHigh","waveguide_model3_flat_long"];
vers = 3;
save_folder = 'test1';
port_i = 2;

matlab_mesh = load(file_list(vers));
Fem_Init(matlab_mesh.no2xyz, matlab_mesh.ed2no_all, matlab_mesh.fa2no_all);
pMtx_ed2no = ProjSol2Nodes_Assemble(matlab_mesh.no2xyz, matlab_mesh.el2no);

max1 = 0; max2 = 0; max3 = 0; max4 = 0;
min1 = 1000; min2 = 1000; min3 = 1000; min4 = 1000;

f_list = (0.25:0.30:1.45)*10^9;
font_size = 11.5;
first = true;
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
%     if kc^2 >= k^2
%         be = 1j*sqrt((k^2) - (kc^2));
%         disp('1j*sqrt((k^2) - (kc^2))')
%     else
%         be = sqrt(((kc^2) - k^2));
%         disp('sqrt(((kc^2) - k^2))')
%     end
    if kc^2 <= k^2
        be = sqrt((k^2) - (kc^2));
        disp('1j*sqrt((k^2) - (kc^2))')
    else
        be = -sqrt((k^2) - (kc^2));
        disp('sqrt(((kc^2) - k^2))')
    end
    
    
%     eta = sqrt(mu/eps);
%     Z = ((k*eta)/be);
%     E0 = sqrt(0.5*Z);
    
    A = (sqrt(2)*pi)/((a^(2/3))*sqrt(b)*sqrt(mu)*sqrt(omega)*sqrt(be));
    %(omega*mu*(a^3)*(A^2)*b)/(
    
    x = 0.05; %the anlytical has the (0,0,0) point at the center inted of the middle
    y = 0.025; %the anlytical has the (0,0,0) point at the center inted of the middle
    z = 0:res:0.4;    
    Ex_ana = 1j*(( 1j*omega*mu*n*pi)/((kc^2)*b)) *A* cos((m*pi*x)/a) * sin((n*pi*y)/b) * exp(-1j*be.*z);
    Ey_ana = 1j*((-1j*omega*mu*m*pi)/((kc^2)*a)) *A* sin((m*pi*x)/a) * cos((n*pi*y)/b) * exp(-1j*be.*z);
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
    
    
    figure(1)
    ax = subplot(3,2,1); hold on
    plot(z,abs(Ey_ana),'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    title('Ey Analytical')
    max1 = max([max1 abs(Ey_ana)]);
    min1 = min([min1 abs(Ey_ana)]);
    ylim([min1,max1])
    xlabel('Z cordinate [m]')
    ylabel('E-Field [V/m]')
    
    ax = subplot(3,2,2); hold on
    plot(CSTz,CST_Eyabs,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    title('CST FEM')
    max2 = max([max2 CST_Eyabs']);
    min2 = min([min2 CST_Eyabs']);
    ylim([min2,max2])
    xlabel('Z cordinate [m]')
    ylabel('E-Field [V/m]')
    
    ax = subplot(3,2,3); hold on
    title('CST FIT')
    plot(CSTz_fit,CST_Eyabs_fit,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    xlabel('Z cordinate [m]')
    ylabel('E-Field [V/m]')
    
    max3 = max([max3 CST_Eyabs_fit']);
    min3 = min([min3 CST_Eyabs_fit']);
    ylim([min3,max3])
    ylabel('E-Field [V/m]')
    
    ax = subplot(3,2,4); hold on
    title('Our FEM')
    plot(Zq,matlab_EyAbs,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    xlabel('Z cordinate [m]')
    max4 = max([max4 matlab_EyAbs]);
    min4 = min([min4 matlab_EyAbs]);
    ylim([min4,max4])
    ylabel('E-Field [V/m]')
    
    lh =legend(ax,'Location','SouthOutside');%'Orientation','Horizontal');
    lh.FontSize = 11.5;
    
    set(gcf,'Position',[100 100 700 460])
    
    figure(2)
    ax = subplot(3,2,1); hold on
    plot(z,angle(Ey_ana),'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    title('Ey Analytical')
    %max1 = max([max1 angle(Ey_ana)]);
    %min1 = min([min1 angle(Ey_ana)]);
    %ylim([min1,max1])
    xlabel('Z cordinate [m]')
    ylabel('E-Field [V/m]')
    
    ax = subplot(3,2,2); hold on
    plot(CSTz,CST_Eyphase,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    title('CST FEM')
    %max2 = max([max2 CST_Eyabs']);
    %min2 = min([min2 CST_Eyabs']);
    %ylim([min2,max2])
    xlabel('Z cordinate [m]')
    ylabel('E-Field [V/m]')
    
    ax = subplot(3,2,3); hold on
    title('CST FIT')
    plot(CSTz_fit,CST_Eyphase_fit,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    xlabel('Z cordinate [m]')
    ylabel('E-Field [V/m]')
    ylabel('E-Field [V/m]')
    
    ax = subplot(3,2,4); hold on
    title('Our FEM')
    plot(Zq,matlab_EyPhase,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    xlabel('Z cordinate [m]')
    ylabel('E-Field [V/m]')
    
    lh =legend(ax,'Location','SouthOutside');%'Orientation','Horizontal');
    lh.FontSize = 11.5;
    
    set(gcf,'Position',[100 100 700 460])
   
end
