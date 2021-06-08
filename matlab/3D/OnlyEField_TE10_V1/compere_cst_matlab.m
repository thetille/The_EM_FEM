clear
close all

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultAxesTickLabelInterpreter', 'none')

file_list = ["cylinder_waveguide2", "waveguide_model3 - simple"...
,"waveguide_model3_flat","mesh_cylinder_R0"...
,"waveguide_model3_flathigh","waveguide_model3_wired"...
,"waveguide_model3_highHigh","waveguide_model3_flat_long"];
vers = 3;
save_folder = 'test1';
port_i = 1;

matlab_mesh = load(file_list(vers));
Fem_Init(matlab_mesh.no2xyz, matlab_mesh.ed2no_all, matlab_mesh.fa2no_all);
pMtx_ed2no = ProjSol2Nodes_Assemble(matlab_mesh.no2xyz, matlab_mesh.el2no);

f_list = (0.1:0.15:1.45)*10^9;
font_size = 11.5;
for f = f_list%0.75*10^9
    if f < 0.75*10^9
        figure_offset = 0;
        cut_off = 'cut';
    else
       figure_offset = 24; 
       cut_off = '';
    end
    %sgtitle(sprintf('f: %f GHz',f*10^-9))
    
    %%%%% Analytical solution %%%%%%
    res = -0.01;
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
    
    
    x = 0.05;
    y = 0.025;
    z = 0:res:-0.4;    
    Ex_ana = (( 1j*omega*mu*n*pi)/((kc^2)*b)) *A* cos((m*pi*x)/a) * sin((n*pi*y)/b) * exp(-1j*be.*z);
    Ey_ana = ((-1j*omega*mu*m*pi)/((kc^2)*a)) *A* sin((m*pi*x)/a) * cos((n*pi*y)/b) * exp(-1j*be.*z);
    Ez_ana = zeros(size(Ey_ana));
    
    
    %%%% plot analytical solution %%%%%
    
    figure(1+figure_offset),hold on
    %subplot(2,3,1), hold on
    plot(z,abs(Ex_ana),'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    title('Magnitude Analytical X')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(2+figure_offset),hold on
    %subplot(2,3,4), hold on
    plot(z,angle(Ex_ana),'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    ylim([-pi, pi])
    set(gca, 'XDir','reverse')
    title('Phase Analytical X')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
    figure(3+figure_offset),hold on
    %subplot(2,3,2), hold on
    plot(z,abs(Ey_ana),'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    title('Magnitude Analytical Y')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    figure(4+figure_offset),hold on
    %subplot(2,3,5), hold on
    plot(z,angle(Ey_ana),'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    ylim([-pi, pi])
    set(gca, 'XDir','reverse')
    title('Phase Analytical Y')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
    figure(5+figure_offset),hold on
    %subplot(2,3,3), hold on
    plot(z,abs(Ez_ana),'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    title('Magnitude Analytical Z')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(6+figure_offset),hold on
    %subplot(2,3,6), hold on
    plot(z,angle(Ez_ana),'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    ylim([-pi, pi])
    set(gca, 'XDir','reverse')
    title('Phase Analytical Z')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
    

    %%%%%% loade CST FEM data %%%%%
    mode = 1;
    if mod(f*10^-9,1) < 0.0001
        if mode == 1
            filename = sprintf('test/cst_res/e-field (f=%.0f) [1].txt',f*10^-9);
        else
            filename = sprintf('test/cst_res_flat/e-field (f=%.0f) [1(1)].txt',f*10^-9);
        end
    elseif mod(f*10^-9,0.1) > 0.01
        if mode == 1
            filename = sprintf('test/cst_res/e-field (f=%.2f) [1].txt',f*10^-9);
        else
            filename = sprintf('test/cst_res_flat/e-field (f=%.2f) [1(1)].txt',f*10^-9);
        end
    else
        if mode == 1
            filename = sprintf('test/cst_res/e-field (f=%.1f) [1].txt',f*10^-9);
        else
            filename = sprintf('test/cst_res_flat/e-field (f=%.1f) [1(1)].txt',f*10^-9);
        end
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

    figure(7+figure_offset),hold on
    %subplot(2,3,1), hold on
    plot(CSTz,CST_Exabs,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    title('Magnitude CST FEM X')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(8+figure_offset),hold on
    %subplot(2,3,4), hold on
    plot(CSTz,CST_Exphase,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    ylim([-pi, pi])
    title('Phase CST FEM X')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
    figure(9+figure_offset),hold on
    %subplot(2,3,2), hold on
    plot(CSTz,CST_Eyabs,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    title('Magnitude CST FEM Y')
    legend()
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(10+figure_offset),hold on
    %subplot(2,3,5), hold on
    plot(CSTz,CST_Eyphase,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    ylim([-pi, pi])
    title('Phase CST FEM Y')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
    figure(11+figure_offset),hold on
    %subplot(2,3,3), hold on
    plot(CSTz,CST_Ezabs,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    title('Magnitude CST FEM Z')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(12+figure_offset),hold on
    %subplot(2,3,6), hold on
    plot(CSTz,CST_Ezphase,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    ylim([-pi, pi])
    title('Phase CST FEM Z')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
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

    figure(13+figure_offset),hold on
    %subplot(2,3,1), hold on
    plot(CSTz_fit,CST_Exabs_fit,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    title('Magnitude CST FIT X')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(14+figure_offset),hold on
    %subplot(2,3,4), hold on
    plot(CSTz_fit,CST_Exphase_fit,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    ylim([-pi, pi])
    title('Phase CST FIT X')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
    figure(15+figure_offset),hold on
    %subplot(2,3,2), hold on
    plot(CSTz_fit,CST_Eyabs_fit,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    title('Magnitude CST FIT Y')
    legend()
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(16+figure_offset),hold on
    %subplot(2,3,5), hold on
    plot(CSTz_fit,CST_Eyphase_fit,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    ylim([-pi, pi])
    title('Phase CST FIT Y')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
    figure(17+figure_offset),hold on
    %subplot(2,3,3), hold on
    plot(CSTz_fit,CST_Ezabs_fit,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    title('Magnitude CST FIT Z')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(18+figure_offset),hold on
    %subplot(2,3,6), hold on
    plot(CSTz_fit,CST_Ezphase_fit,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    set(gca, 'XDir','reverse')
    ylim([-pi, pi])
    title('Phase CST FIT Z')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')
    
    
    %%%%% load matlab data %%%%%
    filename = sprintf('res/%s/%s/E_filds_%d_f_%.0f.mat',file_list(vers),save_folder,port_i,f*10^-6);
    matlab_data = load(filename);
    
    
    %%% interpolate Matlba data %%%%%
    
    %res2 = -res;
    Zq = max(matlab_data.no2xyz(3,:)):res:min(matlab_data.no2xyz(3,:));
    Yq = zeros(size(Zq));
    Xq = zeros(size(Zq));
    
    exFld_all = pMtx_ed2no.xc*matlab_data.eFld_all;
    eyFld_all = pMtx_ed2no.yc*matlab_data.eFld_all;
    ezFld_all = pMtx_ed2no.zc*matlab_data.eFld_all;
    
    matlab_ExRe = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),abs(exFld_all*(-1i)),Xq,Yq,Zq);
    matlab_ExIm = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),angle(exFld_all*(-1i)),Xq,Yq,Zq);
    
    matlab_EyRe = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),abs(eyFld_all*(-1i)),Xq,Yq,Zq);
    matlab_EyIm = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),angle(eyFld_all*(-1i)),Xq,Yq,Zq);
    
    matlab_EzRe = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),abs(ezFld_all*(-1i)),Xq,Yq,Zq);
    matlab_EzIm = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),angle(ezFld_all*(-1i)),Xq,Yq,Zq);
    
    figure(19+figure_offset),hold on
    %subplot(2,3,1), hold on
    plot(Zq,matlab_ExRe,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    title('Magnitude Matlab X')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    figure(20+figure_offset),hold on
    %subplot(2,3,4), hold on
    plot(Zq,matlab_ExIm,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    ylim([-pi, pi])
    title('Phase Matlab X')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
    
    figure(21+figure_offset),hold on
    %subplot(2,3,2), hold on
    plot(Zq,matlab_EyRe,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    title('Magnitude Matlab Y')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(22+figure_offset),hold on
    %subplot(2,3,5), hold on
    plot(Zq,matlab_EyIm,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    ylim([-pi, pi])
    title('Phase Matlab Y')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')

    
    figure(23+figure_offset),hold on
    %subplot(2,3,3), hold on
    plot(Zq,matlab_EzRe,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    title('Magnitude Matlab Z')
    xlabel('Z coordinate [m]')
    ylabel('E-Field [V/m]')

    
    figure(24+figure_offset),hold on
    %subplot(2,3,6), hold on
    plot(Zq,matlab_EzIm,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    ylim([-pi, pi])
    title('Phase Matlab Z')
    xlabel('Z coordinate [m]')
    ylabel('E-Field Phase [rad]')
end


for i = 1:2*figure_offset
    if i < figure_offset, cut_off = 'cut';
    else; cut_off = ''; end
    figure(i)
    lg = legend();
    lg.FontSize = font_size;
    lg.Location = 'eastoutside';
    fil = get(gca,'title');
    fil = get(fil,'string');
    fil = strrep(fil,' ','_');
    fil = sprintf('figures/%s%s',fil,cut_off);
    set(gcf,'Position',[100 100 600 150])
    saveas(gcf,fil,'svg') 
end

close all