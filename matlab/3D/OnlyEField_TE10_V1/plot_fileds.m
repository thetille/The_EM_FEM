%function plot_fileds()%(eFld_all,ed2no_boundery,no2xyz)
clear
file_list = ["cylinder_waveguide2", "waveguide_model3 - simple"...
        ,"waveguide_model3_flat","mesh_cylinder_R0"...
        ,"waveguide_model3_highres","waveguide_model3_wired"...
        ,"waveguide_model3_highHigh","waveguide_model3_flat_long"];
vers = 3;
save_folder = 'test5';
f = 1.2*10^9;

load(sprintf('res/%s/%s/Sparamters',file_list(vers),save_folder))
figure(10), clf
subplot(2,1,1)
plot(f_list*10^(-9),abs(S_par(:,1)),'DisplayName','S11')
hold on
plot(f_list*10^(-9),abs(S_par(:,2)),'DisplayName','S12')
title('abs')
legend()
%ylim([0,2])

a = 0.2;
b = 0.1;
c0 = 299792458; % speed of light in vacuum
w = f_list*2*pi;
k0 = w.*(1/c0);
k_z10 = sqrt((k0.^2)-(pi/a)^2);
Z = (w.*(pi*4*10^(-7))./k_z10);
E0 = sqrt(0.5.*Z);
gamma = 1j*k_z10;

T = ( (2*exp(1j.*k_z10*0.4)) ./ (a*b.*E0) );
s_Analytic = -2.*k_z10.*E0*0.01; 
plot(f_list*10^(-9),abs(s_Analytic*3.5))

subplot(2,1,2)
plot(f_list*10^(-9),angle(S_par(:,1)),'DisplayName','S11')
hold on
plot(f_list*10^(-9),angle(S_par(:,2)),'DisplayName','S12')
title('angle')
legend()

%ylim([0,2])


figure(2), clf


subplot(2,1,1)

plot(f_list*10^(-9),real(S_par(:,1)),'DisplayName','S11')
hold on
plot(f_list*10^(-9),imag(S_par(:,1)),'DisplayName','S11')
legend()
ylim([0,2])
title('real, imag, S11')
subplot(2,1,2)
plot(f_list*10^(-9),real(S_par(:,2)),'DisplayName','S12')
hold on
plot(f_list*10^(-9),imag(S_par(:,2)),'DisplayName','S12')
title('real, imag, S12')
legend()
ylim([0,2])

figure(3), clf;
filename = sprintf('res/%s/%s/E_filds_f_%.0f',file_list(vers),save_folder,f*10^-6);
%filename = 'res/test_save_bMtx';
load(filename);


pMtx_ed2no = ProjSol2Nodes_Assemble(no2xyz, el2no);

exFld_all = pMtx_ed2no.xc*eFld_all;
eyFld_all = pMtx_ed2no.yc*eFld_all;
ezFld_all = pMtx_ed2no.zc*eFld_all;

subplot(1,2,1)

for edIdx = 1:size(ed2no_boundery,2)
    noTmp = ed2no_boundery(:,edIdx);
    xyzTmp = no2xyz(:,noTmp);
    plot3(xyzTmp(1,:), xyzTmp(2,:), xyzTmp(3,:), ...
        'Color', 0.5*[1 1 1])
end

exViz = real(exFld_all(:).');
eyViz = real(eyFld_all(:).');
ezViz = real(ezFld_all(:).');
quiverC3D(no2xyz(1,:)', no2xyz(2,:)', no2xyz(3,:)', ...
    exViz', eyViz', ezViz',1)
set(gca,'Color',[0,27/100,55/100])
axis equal
%axis off
view(40,-20)

subplot(1,2,2), hold on

for edIdx = 1:size(ed2no_boundery,2)
    noTmp = ed2no_boundery(:,edIdx);
    xyzTmp = no2xyz(:,noTmp);
    plot3(xyzTmp(1,:), xyzTmp(2,:), xyzTmp(3,:), ...
        'Color', 0.5*[1 1 1])
end

exViz = imag(exFld_all(:).');
eyViz = imag(eyFld_all(:).');
ezViz = imag(ezFld_all(:).');
quiverC3D(no2xyz(1,:)', no2xyz(2,:)', no2xyz(3,:)', ...
    exViz', eyViz', ezViz',1)
set(gca,'Color',[0,27/100,55/100])

axis equal
%axis off
view(40,-20)

X =  no2xyz(1,:);
Y =  no2xyz(2,:);
Z =  no2xyz(3,:);

V = imag(exFld_all(:))+imag(eyFld_all(:));%+imag(ezFld_all(:));
%V = exFld_all(:).*conj(exFld_all(:))+eyFld_all(:).*conj(eyFld_all(:))+ezFld_all(:).*conj(ezFld_all(:));


res = [0.005,0.005,0.005];
[Xq,Yq,Zq] = meshgrid(-0.1:res(1):0.1, -.1:res(2):.1 ,0:res(3):0.4);
Vq = griddata(X,Y,Z,V,Xq,Yq,Zq);

figure(5), clf;

slice(Xq,Yq,Zq,Vq,[-0.1 0.1],[-0.05 0.05],[0 0.4]);
shading flat
colorbar()
pbaspect([1 1 1])
view(40,-14+(rand(1)*0.1))
title(sprintf('f: %.2f',f*10^-9))
%end