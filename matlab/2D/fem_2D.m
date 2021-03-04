% Physical constants
mu0 = 4*pi*1e-7; % Permeability in vacuum
c0 = 299792456; % Speed of light in vacuum
eps0 = 1/(mu0*c0*c0); % Permittivity in vacuum
% Voltage between inner and outer conductor.
U = 1;
% Read the grid from the file ’unimesh0.mat’.
% This file contains the variables no2xy, el2no, noInt, noExt
load unimesh0
noNum = size(no2xy,2);
elNum = size(el2no,2);
% Scale the domain to measure 2cm x 2cm.
% The initial mesh fitted the unit square:
% -1 < x < 1 and -1 < y < 1.
no2xy = 1e-2*no2xy;
% Assemble the matrix A and vector b.
A = zeros(noNum);
b = zeros(noNum,1);
for elIdx = 1:elNum
    % Get the nodes and their coordinates
    % for the element ’elIdx’.
    no = el2no(:,elIdx);
    xy = no2xy(:,no);
    % Compute the element matrix and add
    % the contribution to the global matrix.
    A_el = CmpElMtx(xy);
    A(no,no) = A(no,no) + A_el;
end
% Get the indices of the nodes.
no_ess = union(noInt, noExt);
no_all = 1:noNum;
no_nat = setdiff(no_all, no_ess);
% Pick out the parts of the matrix and the vectors
% needed to solve the problem.
A_ess = A(no_nat,no_ess);
A_nat = A(no_nat,no_nat);
b = b(no_nat);
z = zeros(length(no_all),1);
z(noInt) = U*ones(length(noInt),1);
z_ess = z(no_ess);
% Solve the system of linear equations.
z_nat = A_nat\(b - A_ess*z_ess);
% Build up the total solution.
z = zeros(length(no_all),1);
z(no_ess) = z_ess;
z(no_nat) = z_nat;
% Compute the capacitance.
W = 0.5*eps0*(z'*A*z);
C=2*W/U^2;
disp(["C per unit length [pF/m] = " num2str(C/1e-12)])
el2ed = unique(el2no);



%%% plot field %%%
plotfield(no2xy, el2no, el2ed, z_nat);