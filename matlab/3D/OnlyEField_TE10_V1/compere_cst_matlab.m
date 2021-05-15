clear

file_list = ["cylinder_waveguide2", "waveguide_model3 - simple"...
,"waveguide_model3_flat","mesh_cylinder_R0"...
,"waveguide_model3_flathigh","waveguide_model3_wired"...
,"waveguide_model3_highHigh","waveguide_model3_flat_long"];
vers = 3;
save_folder = 'test5';

%addpath("C:\Users\benja\Sync\AAU\10. semester\github\The_EM_FEM\matlab\3D\OnlyEField_TE10_V1")
matlab_mesh = load(file_list(vers));
Fem_Init(matlab_mesh.no2xyz, matlab_mesh.ed2no_all, matlab_mesh.fa2no_all);
pMtx_ed2no = ProjSol2Nodes_Assemble(matlab_mesh.no2xyz, matlab_mesh.el2no);

f_list = (0.6:0.05:1.2)*10^9;
figure(1), clf
for f = f_list%0.75*10^9
    sgtitle(sprintf('f: %f GHz',f*10^-9))

    %%%%%% loade CST data %%%%%
    mode = 3;
    if mod(f*10^-9,1) < 0.0001
        if mode == 1
            filename = sprintf('res/cst_res2/e-field (f=%.0f) [].txt',f*10^-9);
        else
            filename = sprintf('test/cst_res_flat/e-field (f=%.0f) [1(1)].txt',f*10^-9);
        end
    elseif mod(f*10^-9,0.1) > 0.01
        if mode == 1
            filename = sprintf('test/cst_res2/e-field (f=%.2f) [1].txt',f*10^-9);
        else
            filename = sprintf('test/cst_res_flat/e-field (f=%.2f) [1(1)].txt',f*10^-9);
        end
    else
        if mode == 1
            filename = sprintf('test/cst_res2/e-field (f=%.1f) [1].txt',f*10^-9);
        else
            filename = sprintf('test/cst_res_flat/e-field (f=%.1f) [1(1)].txt',f*10^-9);
        end
    end
    delimiterIn = ' ';
    headerlinesIn = 2;
    CST_data = importdata(filename,delimiterIn,headerlinesIn);
    CST_data = CST_data.data; %x,y,z, ExRe,EyRe,EzRe, ExIm,EyIm,EzIm 
    
    %%%%% load matlab data %%%%%
    filename = sprintf('res/%s/%s/E_filds_f_%.0f.mat',file_list(vers),save_folder,f*10^-6);
    matlab_data = load(filename);
    
    
    %%%% interpolate and get center line for CST data %%%%%
    res = 0.001;
    Zq = min(CST_data(:,3)):res:max(CST_data(:,3));
    Yq = zeros(size(Zq));
    Xq = zeros(size(Zq));
    
    CST_ExRe = griddata(CST_data(:,1),CST_data(:,2),CST_data(:,3),abs(CST_data(:,4)+CST_data(:,4+3)*1i),Xq,Yq,Zq);
    CST_ExIm = griddata(CST_data(:,1),CST_data(:,2),CST_data(:,3),angle(CST_data(:,4)+CST_data(:,4+3)*1i),Xq,Yq,Zq);
    
    CST_EyRe = griddata(CST_data(:,1),CST_data(:,2),CST_data(:,3),abs(CST_data(:,5)+CST_data(:,5+3)*1i),Xq,Yq,Zq);
    CST_EyIm = griddata(CST_data(:,1),CST_data(:,2),CST_data(:,3),angle(CST_data(:,5)+CST_data(:,5+3)*1i),Xq,Yq,Zq);
    
    CST_EzRe = griddata(CST_data(:,1),CST_data(:,2),CST_data(:,3),abs(CST_data(:,6)+CST_data(:,6+3)*1i),Xq,Yq,Zq);
    CST_EzIm = griddata(CST_data(:,1),CST_data(:,2),CST_data(:,3),angle(CST_data(:,6)+CST_data(:,6+3)*1i),Xq,Yq,Zq);
    
    subplot(4,3,1), hold on
    plot(Zq,CST_ExRe)
    title('abs cst x')
    subplot(4,3,4), hold on
    plot(Zq,CST_ExIm)
    title('angle cst x')
    
    subplot(4,3,2), hold on
    plot(Zq,CST_EyRe,'DisplayName',sprintf('f: %.2f GHz',f*10^-9))
    title('abs cst y')
    
    subplot(4,3,5), hold on
    plot(Zq,CST_EyIm)
    title('angle cst y')
    
    subplot(4,3,3), hold on
    plot(Zq,CST_EzRe)
    title('abs cst z')
    subplot(4,3,6), hold on
    plot(Zq,CST_EzIm)
    title('angle cst z')
    
    %%% interpolate Matlba data %%%%%

  
    
    
    
    Zq = min(matlab_data.no2xyz(3,:)):res:max(matlab_data.no2xyz(3,:));
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
    
    subplot(4,3,7), hold on
    plot(Zq,matlab_ExRe)
    title('abs matlab x')
    subplot(4,3,10), hold on
    plot(Zq,matlab_ExIm)
    title('angle matlab x')
    
    subplot(4,3,8), hold on
    plot(Zq,matlab_EyRe)
    title('abs matlab y')
    subplot(4,3,11), hold on
    plot(Zq,matlab_EyIm)
    title('angle matlab y')
    
    subplot(4,3,9), hold on
    plot(Zq,matlab_EzRe)
    title('abs matlab z')
    subplot(4,3,12), hold on
    plot(Zq,matlab_EzIm)
    title('angle matlab z')
    
%     figure(2), clf
%     exViz = real(exFld_all(:).');
%     eyViz = real(eyFld_all(:).');
%     ezViz = real(ezFld_all(:).');
%     quiverC3D(matlab_data.no2xyz(1,:)', matlab_data.no2xyz(2,:)', matlab_data.no2xyz(3,:)', ...
%         exViz', eyViz', ezViz',1)
%     set(gca,'Color',[0,27/100,55/100])
    
end
subplot(4,3,2)
legend();