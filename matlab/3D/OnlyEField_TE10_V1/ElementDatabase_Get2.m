function idx_DB = ElementDatabase_Get(name_ipt, data_ipt)

% Call:
%    ElementDatabase_Init(name_ipt, data_ipt)
% Arguments:
%    name_ipt = name of the database
%    data_ipt = data to store in the database [length(dim) x # entries]
% Returns:
%    idx_DB = indices in the database that correspond to data_ipt
% Comments:
%    Can be used to store edges and faces. Given that the number of data
%    entries is referred to as daNum, we have the following examples
%
%    1) Edges: size(data_ipt) = [2 daNum]
%    2) Triangular faces: size(data_ipt) = [3 daNum]
%    3) Quadrilateral faces: size(data_ipt) = [4 daNum]

% Access the global variable
str = sprintf('global gbl___%s', name_ipt);
eval(str)

% Get the dimensions and sizes of the database
str = sprintf('size_DB = gbl___%s.size;', name_ipt);
eval(str)

% Get the hash-table keys of the database
str = sprintf('htk_DB = gbl___%s.htk;', name_ipt);
eval(str)

% Check dimensions
if length(size_DB) ~= size(data_ipt,1)
    error('Incorrect dimension for the data')
end

% Check data
for dimIdx = 1:length(size_DB)
    minTmp = min(data_ipt(dimIdx));
    maxTmp = max(data_ipt(dimIdx));
    if (minTmp < 0 || maxTmp > size_DB(dimIdx))
        error('Data out of bounds')
    end
end

% Compute the hash-table key from the data
data_ipt = sort(data_ipt);
htk_ipt = data_ipt(1,:);
jmpTmp = size_DB(1);
for dimIdx = 2:length(size_DB)
    htk_ipt = htk_ipt + jmpTmp*(data_ipt(dimIdx,:)-1);
    jmpTmp = jmpTmp*size_DB(dimIdx);
end

[~, idx_DB, idx_ipt] = intersect(htk_DB, htk_ipt);
[~, idx_sort] = sort(idx_ipt);

idx_ipt = idx_ipt(idx_sort);
idx_DB = idx_DB(idx_sort)';

% ------------------------------------------------------------------
% Sanity check (Adds about 10% execution time...)
% ------------------------------------------------------------------
return

% Get the hash-table keys of the database
str = sprintf('data_DB = gbl___%s.data;', name_ipt);
eval(str)

if norm(data_DB(:,idx_DB) - data_ipt) > 0
    keyboard
end

num = ElementDatabase_Cardinal(name_ipt);
if max(data_DB) > num
    keyboard
end