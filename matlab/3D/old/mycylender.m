R = 0.125 % radius of cylinder
H = 0.4 % hight of cylinder
model = createpde;
gm = multicylinder(R,H);
model.Geometry = gm;
mesh_default = generateMesh(model, 'GeometricOrder', 'linear');
figure
pdemesh(mesh_default)

no2xyz = mesh_default.Nodes; % coordinates of the nodes
el2no = mesh_default.Elements; % nodes of the tetrahedrons
el2ma = 2.5 * ones(1,size(mesh_default.Elements,2)); % material indices of the tetrahedrons
ed2no_pec = [] ; % nodes of the edges that are located on the surface S (PEC)
for i = 1:size(no2xyz,2)
    if (no2xyz(1,i)^2+no2xyz(2,i)^2) > R^2-0.01 || no2xyz(3,i)^2 > H^2-0.1
       ed2no_pec(:,i) = no2xyz(:,i); 
    end
end
for i = 1:size(ed2no_pec,2)
    if ed2no_pec(2,i) < 0
        scatter3(ed2no_pec(1,i),ed2no_pec(2,i),ed2no_pec(3,i));
    end
end
%ed2no_all = ; % nodes of all edges in the mesh
%fa2no_all = ; % nodes of all faces (i.e. triangles) in the mesh

