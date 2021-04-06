function crd_DB = ElementDatabase_Cardinal(name_ipt)

% Call:
%    ElementDatabase_Cardinal(name_ipt)
% Arguments:
%    name_ipt = name of the database
% Returns:
%    crd_DB = cardinal (number of entries)
% Comments:
%    -

% Access the global variable
str = sprintf('global gbl___%s', name_ipt);
eval(str)

% Get the hash-table keys of the database
str = sprintf('htk_DB = gbl___%s.htk;', name_ipt);
eval(str)

% Compute the cardinal of the database
crd_DB = size(htk_DB,2);