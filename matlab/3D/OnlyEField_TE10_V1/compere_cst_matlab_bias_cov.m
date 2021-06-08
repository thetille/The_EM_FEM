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

f_list = (0.25:0.30:1.45)*10^9;
font_size = 11.5;

% Eyabs_ana_mat = [];
% Eyabs_ana_CSTFEM = [];
% Eyabs_ana_CSTFit = [];
% Eyabs_mat_CST = [];
% Eyabs_mat_CSTFit = [];
% Eyabs_CSTFEM_CSTFit = [];

Eyabs_ana_CSTFIT = [];
Eyabs_ana_CSTFEM = [];
Eyabs_ana_mat = [];

for f = f_list
    if f < 0.75*10^9
        figure_offset = 0;
        cut_off = 'cut';
    else
       figure_offset = 18; 
       cut_off = '';
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
    be = -sqrt((k^2) - (kc^2));
    
%     eta = sqrt(mu/eps);
%     Z = ((k*eta)/be);
%     E0 = sqrt(0.5*Z);
    
    
    A = (sqrt(2)*pi)/((a^(2/3))*sqrt(b)*sqrt(mu)*sqrt(omega)*sqrt(be));
    
    %(omega*mu*(a^3)*(A^2)*b)/(
    
    x = 0.05;
    y = 0.025;
    z = 0:res:0.4;    
    Ex_ana = (( 1j*omega*mu*n*pi)/((kc^2)*b)) *A* cos((m*pi*x)/a) * sin((n*pi*y)/b) * exp(-1j*be.*z);
    Ey_ana = ((-1j*omega*mu*m*pi)/((kc^2)*a)) *A* sin((m*pi*x)/a) * cos((n*pi*y)/b) * exp(-1j*be.*z);
    Ez_ana = zeros(size(Ey_ana));

    %%%%%% loade CST data %%%%%
    
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
      
    %%%% CST data %%%

    CST_Exabs = abs(CST_data(:,4)+CST_data(:,4+3)*1i);
    CST_Exphase = angle(CST_data(:,4)+CST_data(:,4+3)*1i);
    
    CST_Eyabs = abs(CST_data(:,5)+CST_data(:,5+3)*1i);
    CST_Eyphase = angle(CST_data(:,5)+CST_data(:,5+3)*1i);
    
    CST_Ezabs = abs(CST_data(:,6)+CST_data(:,6+3)*1i);
    CST_Ezphase = angle(CST_data(:,6)+CST_data(:,6+3)*1i);

    CSTz = CST_data(:,3);
    
    
    
    %%%% loade CST FIT data %%%%%
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
      
    %%%% CST data %%%

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
    
    
%     Eyabs_ana_mat =         [Eyabs_ana_mat, mean(abs(Ey_ana)./matlab_EyAbs)];
%     Eyabs_ana_CSTFEM =      [Eyabs_ana_CSTFEM, mean(abs(Ey_ana)./CST_Eyabs')];
%     Eyabs_ana_CSTFit =      [Eyabs_ana_CSTFit, mean(abs(Ey_ana)./CST_Eyabs_fit')];
%     Eyabs_mat_CST =         [Eyabs_mat_CST, mean(matlab_EyAbs./CST_Eyabs')];
%     Eyabs_mat_CSTFit =      [Eyabs_mat_CSTFit, mean(matlab_EyAbs./CST_Eyabs_fit')];
%     Eyabs_CSTFEM_CSTFit =   [Eyabs_CSTFEM_CSTFit, mean(CST_Eyabs'./CST_Eyabs_fit')];
    

    Eyabs_ana_CSTFIT =  [Eyabs_ana_CSTFIT, mean(abs(Ey_ana)./CST_Eyabs_fit')];
    Eyabs_ana_CSTFEM =  [Eyabs_ana_CSTFEM, mean(abs(Ey_ana)./CST_Eyabs')];
    Eyabs_ana_mat =     [Eyabs_ana_mat,    mean(abs(Ey_ana)./matlab_EyAbs)];
end

font_size = 11.5;

save('normelization.mat','Eyabs_ana_CSTFIT','Eyabs_ana_CSTFEM','Eyabs_ana_mat')

% figure(1), hold on
% subplot(2,3,1)
% plot(f_list*10^-9,abs(Eyabs_ana_mat))
% title('Analytical Vs. Our FEM')
% subplot(2,3,2)
% plot(f_list*10^-9,abs(Eyabs_ana_CSTFEM))
% title('Analytical Vs. CST FEM')
% subplot(2,3,3)
% plot(f_list*10^-9,abs(Eyabs_ana_CSTFit))
% title('Analytical Vs. CST FIT')
% subplot(2,3,4)
% plot(f_list*10^-9,abs(Eyabs_mat_CST))
% title('Our FEM Vs. CST FEM')
% subplot(2,3,5)
% plot(f_list*10^-9,abs(Eyabs_mat_CSTFit))
% title('Our FEM Vs. CST FIT')
% subplot(2,3,6)
% plot(f_list*10^-9,abs(Eyabs_CSTFEM_CSTFit))
% title('CST FEM Vs. CST FIT')

figure(1)
subplot(1,3,1), hold on
stem(f_list*10^-9,(Eyabs_ana_CSTFIT),'color',[0,0,0])
ylim([0.24,0.27])
title('Analytical Vs. CST FIT')
subplot(1,3,2), hold on
stem(f_list*10^-9,(Eyabs_ana_CSTFEM),'color',[0,0,0])
ylim([0.262,0.267])
title('Analytical Vs. CST FEM')
subplot(1,3,3), hold on
stem(f_list*10^-9,(Eyabs_ana_mat),'color',[0,0,0])
ylim([3.68,3.74])
title('Analytical Vs. Our FEM')



%title('Absolute difference of X Component')

% figure(2), hold on
% plot(f_list*10^-9,abs(Eyabs_ana_mat),'DisplayName','Analytical Vs. our FEM')
% plot(f_list*10^-9,abs(Eyabs_CST_mat),'DisplayName','CST FEM Vs. our FEM')
% plot(f_list*10^-9,abs(Eyabs_ana_cst),'DisplayName','CST FEM Vs. Analytical')
% plot([0.74948, 0.74948],get(gca,'YLim'),'--','color',[0.7,0.7,0.7],'DisplayName','Cut off frequency')
% title('Absolute difference of Y Component')
% 
% figure(3), hold on
% plot(f_list*10^-9,abs(Ezabs_ana_mat),'DisplayName','Analytical Vs. our FEM')
% plot(f_list*10^-9,abs(Ezabs_CST_mat),'DisplayName','CST FEM Vs. our FEM')
% plot(f_list*10^-9,abs(Ezabs_ana_cst),'DisplayName','CST FEM VS. Analytical')
% 
% plot([0.74948, 0.74948],get(gca,'YLim'),'--','color',[0.7,0.7,0.7],'DisplayName','Cut off frequency')
% title('Absolute difference of Z Component')

figure(1)
for i = 1:3
    subplot(1,3,i)
    xlabel('Freqency GHZ')
    ylabel('Analytical/CST')
    ax = gca;
    ylim(get(ax,'YLim'))%+[-0.000001 0.000001])
    xlim(get(ax,'XLim'))%+[-0.000001 0.000001])
    plot([0.74948, 0.74948],get(ax,'YLim'),'--','color',[0.7,0.7,0.7])
end


% for i = 1:3
%     figure(i)
%     xlabel('Freqency GHZ')
%     ylabel('Difference [V/m]')
%     grid on
%     %set(gca, 'YScale', 'log')
%     lg = legend();
%     lg.FontSize = font_size;
%     lg.Location = 'eastoutside';
%     set(gcf,'Position',[100 100 600 150])
%     fil = get(gca,'title');
%     fil = get(fil,'string');
%     fil = strrep(fil,' ','_');
%     fil = sprintf('figures2/%s',fil);
%     set(gcf,'Position',[100 100 600 150])
%     saveas(gcf,fil,'svg')
% end


% figure(4), hold on
% plot(f_list*10^-9,Exabs_ana_mat_cov,'DisplayName','Ex analytical Vs. FEM implementation')
% plot(f_list*10^-9,Exabs_CST_mat_cov,'DisplayName','Ex CST Vs. FEM implementation')
% plot(f_list*10^-9,Exabs_ana_cst_cov,'DisplayName','Ex Analytical Vs. CST')
% plot([0.74948, 0.74948],get(gca,'YLim'),'--','color',[0.7,0.7,0.7],'DisplayName','Cut off frequency')
% title('Covariance Comparison of X Component')
% 
% figure(5), hold on
% plot(f_list*10^-9,Eyabs_ana_mat_cov,'DisplayName','Ey analytical Vs. FEM implementation')
% plot(f_list*10^-9,Eyabs_CST_mat_cov,'DisplayName','Ey CST Vs. FEM implementation')
% plot(f_list*10^-9,Eyabs_ana_cst_cov,'DisplayName','Ey Analytical Vs. CST')
% plot([0.74948, 0.74948],get(gca,'YLim'),'--','color',[0.7,0.7,0.7],'DisplayName','Cut off frequency')
% title('Covariance Comparison of Y Component')
% 
% figure(6), hold on
% plot(f_list*10^-9,Ezabs_ana_mat_cov,'DisplayName','Ez analytical Vs. FEM implementation')
% plot(f_list*10^-9,Ezabs_CST_mat_cov,'DisplayName','Ez CST Vs. FEM implementation')
% plot(f_list*10^-9,Ezabs_ana_cst_cov,'DisplayName','Ez Analytical Vs. CST')
% plot([0.74948, 0.74948],get(gca,'YLim'),'--','color',[0.7,0.7,0.7],'DisplayName','Cut off frequency')
% title('Covariance Comparison of Z Component')
% 
% for i = 1:3
%     figure(i+3)
%     xlabel('Freqency GHZ')
%     ylabel('Covariance')
%     lg = legend();
%     lg.FontSize = font_size;
%     lg.Location = 'eastoutside';
%     set(gcf,'Position',[100 100 600 150])
%     fil = get(gca,'title');
%     fil = get(fil,'string');
%     fil = strrep(fil,' ','_');
%     fil = sprintf('figures2/%s',fil);
%     set(gcf,'Position',[100 100 600 150])
%     saveas(gcf,fil,'svg')  
% end


% figure(7), hold on
% plot(f_list*10^-9,Eyabs_ana_mat_cof,'DisplayName','Ey analytical Vs. FEM implementation')
% plot(f_list*10^-9,Eyabs_CST_mat_cof,'DisplayName','Ey CST Vs. FEM implementation')
% plot(f_list*10^-9,Eyabs_ana_cst_cof,'DisplayName','Ey Analytical Vs. CST')
% plot([0.74948, 0.74948],get(gca,'YLim'),'--','color',[0.7,0.7,0.7],'DisplayName','Cut off frequency')
% title('Correlation Comparison of X Component')
% 
% figure(8), hold on
% plot(f_list*10^-9,Eyabs_ana_mat_cof,'DisplayName','Ey analytical Vs. FEM implementation')
% plot(f_list*10^-9,Eyabs_CST_mat_cof,'DisplayName','Ey CST Vs. FEM implementation')
% plot(f_list*10^-9,Eyabs_ana_cst_cof,'DisplayName','Ey Analytical Vs. CST')
% plot([0.74948, 0.74948],get(gca,'YLim'),'--','color',[0.7,0.7,0.7],'DisplayName','Cut off frequency')
% title('Correlation Comparison of Y Component')
% 
% figure(9), hold on
% plot(f_list*10^-9,Eyabs_ana_mat_cof,'DisplayName','Ey analytical Vs. FEM implementation')
% plot(f_list*10^-9,Eyabs_CST_mat_cof,'DisplayName','Ey CST Vs. FEM implementation')
% plot(f_list*10^-9,Eyabs_ana_cst_cof,'DisplayName','Ey Analytical Vs. CST')
% plot([0.74948, 0.74948],get(gca,'YLim'),'--','color',[0.7,0.7,0.7],'DisplayName','Cut off frequency')
% title('Correlation Comparison of Z Component')
% 
% for i = 1:3
%     figure(i+6)
%     xlabel('Freqency GHZ')
%     ylabel('Correlation')
%     lg = legend();
%     lg.FontSize = font_size;
%     lg.Location = 'eastoutside';
%     set(gcf,'Position',[100 100 600 150])
%     fil = get(gca,'title');
%     fil = get(fil,'string');
%     fil = strrep(fil,' ','_');
%     fil = sprintf('figures2/%s',fil);
%     set(gcf,'Position',[100 100 600 150])
%     saveas(gcf,fil,'svg')  
% end
% 
% function diff = NMSD(x,y)
%     normelized_x =  
% 
% 
% end