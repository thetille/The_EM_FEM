clear


%%% loade data %%%
mesh = 0:10;
f = (0.25:0.30:1.50)*10^9;

addpath('C:\Users\benja\Sync\AAU\10. semester\github\The_EM_FEM\matlab\3D\OnlyEField_TE10_V1')

for mesh_i = 1:length(mesh)
    mesh(mesh_i)
    %loade mesh
    matlab_mesh = load(sprintf('waveguide_model3_flat%d.mat',mesh(mesh_i)));
    Fem_Init(matlab_mesh.no2xyz, matlab_mesh.ed2no_all, matlab_mesh.fa2no_all);
    pMtx_ed2no = ProjSol2Nodes_Assemble(matlab_mesh.no2xyz, matlab_mesh.el2no);

    %loade simulation data
    for fi = 1:length(f)
        filename = sprintf('waveguide_model3_flat%d/test1/E_filds_%d_f_%.0f.mat',mesh(mesh_i),1,f(fi)*10^-6);
        matlab_data = load(filename);
    
        %find the E fileds at nodes
        exFld_all = pMtx_ed2no.xc*matlab_data.eFld_all;
        eyFld_all = pMtx_ed2no.yc*matlab_data.eFld_all;
        ezFld_all = pMtx_ed2no.zc*matlab_data.eFld_all;
    
        %find the E fileds at line
        res = 0.01;
        Zq = min(matlab_data.no2xyz(3,:)):res:max(matlab_data.no2xyz(3,:));
        Yq = ones(size(Zq))*0.025;
        Xq = ones(size(Zq))*-0.05;
    
        matlab_ExAbs = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),abs(exFld_all*(-1i)),Xq,Yq,Zq);
        matlab_ExPhase = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),angle(exFld_all*(-1i)),Xq,Yq,Zq);

        matlab_EyAbs = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),abs(eyFld_all*(-1i)),Xq,Yq,Zq);
        matlab_EyPhase = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),angle(eyFld_all*(-1i)),Xq,Yq,Zq);

        matlab_EzAbs = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),abs(ezFld_all*(-1i)),Xq,Yq,Zq);
        matlab_EzPhase = griddata(matlab_data.no2xyz(1,:),matlab_data.no2xyz(2,:),matlab_data.no2xyz(3,:),angle(ezFld_all*(-1i)),Xq,Yq,Zq);

        matlab_ExAbs = matlab_ExAbs/mean(matlab_ExAbs);
        matlab_EyAbs = matlab_EyAbs/mean(matlab_EyAbs);
        matlab_ExAbs = matlab_EzAbs/mean(matlab_EzAbs);
        
        if mesh_i ~= 1

            differencex(mesh_i,fi) = mean((matlab_ExAbs-matlab_ExAbs_old(fi,:)).^2);
            differencey(mesh_i,fi) = mean((matlab_EyAbs-matlab_EyAbs_old(fi,:)).^2);
            differencez(mesh_i,fi) = mean((matlab_EzAbs-matlab_EzAbs_old(fi,:)).^2);

        end
    
    
        matlab_ExAbs_old(fi,:) = matlab_ExAbs;
        matlab_EyAbs_old(fi,:) = matlab_EyAbs;
        matlab_EzAbs_old(fi,:) = matlab_EzAbs;
    end
end

save('results_mesh_desety')
mesh_tetras = [733, 1714, 2770, 4982, 11639, 19785, 37749, 46831, 61127, 88456, 134590];
figure(1), clf, hold on
for fi = 1:length(f)
    plot(mesh_tetras(mesh(2:end)+1),differencex(2:end,fi),'-','Marker','*','DisplayName',sprintf('f: %.2f',f(fi)*10^-9))
end
xlim([1700,mesh_tetras(mesh_i)])
ylabel('Change in the E-Field Amplitude')
xlabel('Number of Tetrahedrons')
title('Convergence Test for the E-fields X Component')
%xticks(mesh_tetras)
%xticklabels(mesh_tetras)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on
lg = legend();
lg.FontSize = 11.5;

figure(2), clf, hold on
for fi = 1:length(f)
    plot(mesh_tetras(mesh(2:end)+1),differencey(2:end,fi),'Marker','*','DisplayName',sprintf('f: %.2f',f(fi)*10^-9))
end
xlim([1000,mesh_tetras(mesh_i)])
ylabel('Change in the E-Field Amplitude')
xlabel('Number of Tetrahedrons')
title('Convergence Test for the E-fields Y Component')
%xticks(mesh_tetras)
%xticklabels(mesh_tetras)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on
lg = legend();
lg.FontSize = 11.5;

figure(3), clf, hold on
for fi = 1:length(f)
    plot(mesh_tetras(mesh(2:end)+1),differencez(2:end,fi),'Marker','*','DisplayName',sprintf('f: %.2f',f(fi)*10^-9))
end
xlim([1700,mesh_tetras(mesh_i)])
ylabel('Change in the E-Field Amplitude')
xlabel('Number of Tetrahedrons')
title('Convergence Test for the E-fields Z Component')
%xticks(mesh_tetras)
%xticklabels(mesh_tetras)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on
lg = legend();
lg.FontSize = 11.5;
