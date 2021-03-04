% --------------------------------------------------------------
% Plot a 2D vector field described by edge elements
% --------------------------------------------------------------
function plotfield(no2xy, el2no, el2ed, sol)
% Arguments:
% no2xy = x- and y-coordinates of the nodes
% el2no = node indices for all triangles
% el2ed = edge indices for all elements
% sol = Coefficient vector (each entry in the vector
% corresponds to one edge in the mesh)
% Returns:
% Sort the nodes of each element
el2no = sort(el2no);
% Local coordinates for subgrid plotting
phi_1 = [4 3 2 1 0 3 2 1 0 2 1 0 1 0 0]' / 4;
phi_2 = [0 1 2 3 4 0 1 2 3 0 1 2 0 1 0]' / 4;
phi_3 = [0 0 0 0 0 1 1 1 1 2 2 2 3 3 4]' / 4;
% Gradients of the simplex functions
% (constant within each element)
edge1 = no2xy(:,el2no(2,:)) - no2xy(:,el2no(1,:));
edge2 = no2xy(:,el2no(3,:)) - no2xy(:,el2no(1,:));
detJ = edge1(1,:).*edge2(2,:) - edge1(2,:).*edge2(1,:);
grad_phi_2x = edge2(2,:)./ detJ;
grad_phi_2y = -edge2(1,:)./ detJ;
grad_phi_3x = -edge1(2,:)./ detJ;
grad_phi_3y = edge1(1,:)./ detJ;
grad_phi_1x=0- grad_phi_2x - grad_phi_3x;
grad_phi_1y=0- grad_phi_2y - grad_phi_3y;
% Solution values associated to the 1st, 2nd, and
% 3rd edges in each element
sol1 = sol(el2ed(1,:)).';
sol2 = sol(el2ed(2,:)).';
sol3 = sol(el2ed(3,:)).';
% Field values
Ex = phi_1 * ( grad_phi_2x.*sol1 + grad_phi_3x.*sol2) + ...
    phi_2 * (-grad_phi_1x.*sol1 + grad_phi_3x.*sol3) + ...
    phi_3 * (-grad_phi_1x.*sol2 - grad_phi_2x.*sol3);
Ey = phi_1 * ( grad_phi_2y.*sol1 + grad_phi_3y.*sol2) + ...
    phi_2 * (-grad_phi_1y.*sol1 + grad_phi_3y.*sol3) + ...
    phi_3 * (-grad_phi_1y.*sol2 - grad_phi_2y.*sol3);
Hz = (sol1 - sol2 + sol3)./detJ;
% Create subgrid
p1 = no2xy(:,el2no(1,:));
p2 = no2xy(:,el2no(2,:));
p3 = no2xy(:,el2no(3,:));
psub = kron(p1,phi_1') + kron(p2,phi_2') + kron(p3,phi_3');
% Initiate plotting
ih = ishold;
ax = newplot;
% Plot the curl of the field (constant within each element)
patch('faces',el2no','vertices',no2xy','facevertexcdata',Hz(:), ...
    'facecolor',get(ax,'defaultsurfacefacecolor'), ...
    'edgecolor',get(ax,'defaultsurfaceedgecolor'));
axis equal, hold on
% Plot the field itself as arrows
quiver(psub(1,:),psub(2,:),Ex(:)',Ey(:)','k');
% Plot the mesh
xy1 = no2xy(:,el2no(1,:));
xy2 = no2xy(:,el2no(2,:));
xy3 = no2xy(:,el2no(3,:));
xy = [xy1; xy2; xy3; xy1; NaN*xy1];
plot(xy(1:2:end),xy(2:2:end),'k')
% Create a new colormap
mrz = max(abs(Hz(:)));
caxis([-mrz, mrz]);
c = (0:64)'/64; d = [c c ones(size(c))];
colormap([d ;1 1 1; d(end:-1:1,end:-1:1)]);
if ~ih, hold off, end