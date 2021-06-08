%%%%
% The code in this file is from the book 
% Computational electromagnetics by
% Thomas Rylander, Anders Bondeson, Par Ingelstrom
%%%%


function ElementDatabase_Set(name_ipt, data_ipt)

% Call:
%    ElementDatabase_Init(name_ipt, data_ipt)
% Arguments:
%    name_ipt = name of the database
%    data_ipt = data to store in the database [length(dim) x # entries]
% Returns:
%    -
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

% Get the dimension and size of the database
str = sprintf('size_DB = gbl___%s.size;', name_ipt);
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
htkVal = sort(data_ipt);
htkIdx = htkVal(1,:);
jmpTmp = size_DB(1);
for dimIdx = 2:length(size_DB)    
    htkIdx = htkIdx + jmpTmp*(htkVal(dimIdx,:)-1);
    jmpTmp = jmpTmp*size_DB(dimIdx);
end

% Set the data (= hash-table keys)
str = sprintf('gbl___%s.data = htkVal;', name_ipt);
eval(str)

% Set the indices (= hash-table indices in linear dimension)
str = sprintf('gbl___%s.htk = htkIdx;', name_ipt);
eval(str)
