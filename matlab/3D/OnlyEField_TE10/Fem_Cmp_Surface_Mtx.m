% --------------------------------------------------------------
% Compute element matrices for the port surfaces by means of
% numerical integration on the reference element
% --------------------------------------------------------------
function [KElMtx_EE] = ...
    Fem_Cmp_Surface_Mtx(xyz, ma2er, ma2si, k0, ed2no_port1, ed2no_port2) % need to change from hardcode to adaptive
% Argument:
%   xyz = the coordinates of the nodes of the element
%   ma2er = material to permittivity
%   ma2si = material to conductivity
% Return:
%
