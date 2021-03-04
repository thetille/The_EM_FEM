% --------------------------------------------------------------
% Compute element matrix for a triangle and its node basis
% --------------------------------------------------------------
function Ae = CmpElMtx(xy)
% Arguments:
% xy = the coordinates of the nodes of the triangle
% Returns:
% Ae = element matrix corresponding to the Laplace operator
% Edges
s1 = xy(:,3)-xy(:,2);
s2 = xy(:,1)-xy(:,3);
s3 = xy(:,2)-xy(:,1);
% Area of the triangle
Atot = 0.5*(s2(1)*s3(2)-s2(2)*s3(1));
% Check whether area is negative (nodes given counterclockwise)
if (Atot < 0)
error("The nodes of the element given in wrong order")
end
% Compute the gradient of the vectors.
grad_phi1e = [-s1(2);s1(1)]/(2*Atot);
grad_phi2e = [-s2(2);s2(1)]/(2*Atot);
grad_phi3e = [-s3(2);s3(1)]/(2*Atot);
grad_phi = [grad_phi1e grad_phi2e grad_phi3e];
% Compute all the integrals for this particular element.
for iIdx = 1:3
for jIdx = 1:3
Ae(iIdx,jIdx) = grad_phi(:,iIdx)' * grad_phi(:,jIdx) * Atot;
end
end