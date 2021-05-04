%function plot_fileds()%(eFld_all,ed2no_boundery,no2xyz)

file_list = ["cylinder_waveguide2", "waveguide_model3 - simple"...
        ,"waveguide_model3","mesh_cylinder_R0"...
        ,"waveguide_model3_highres","waveguide_model3_wired"...
        ,"waveguide_model3_highHigh"];
vers = 5;

f = 1.1*10^9;
filename = sprintf('res/%s/E_filds_f_%.0f',file_list(vers),f*10^-6);
load(filename);

pMtx_ed2no = ProjSol2Nodes_Assemble(no2xyz, el2no);

exFld_all = pMtx_ed2no.xc*eFld_all;
eyFld_all = pMtx_ed2no.yc*eFld_all;
ezFld_all = pMtx_ed2no.zc*eFld_all;


figure(3), clf;
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

V = abs(exViz)+abs(eyViz)+abs(ezViz);

res = [0.005,0.005,0.005];
[Xq,Yq,Zq] = meshgrid(-0.1:res(1):0.1, -.1:res(2):.1 ,0:res(3):0.4);
Vq = griddata(X,Y,Z,V,Xq,Yq,Zq);

figure(5), clf;
slice(Xq,Yq,Zq,Vq,[-0.1 0.1],[-0.1 0.1],[0 0.4]);
shading flat
pbaspect([1 1 1])
view(40,-14+(rand(1)*0.1))
%end