%%%%
% The code in this file is from the book 
% Computational electromagnetics by
% Thomas Rylander, Anders Bondeson, Par Ingelstrom
%%%%


function ElementDatabase_Init(name_ipt, size_ipt)

% Call:
%    ElementDatabase_Init(name_ipt, size_ipt)
% Arguments:
%    name_ipt = name of the database
%    size_ipt = number of dimensions and their sizes
% Returns:
%    -
% Comments:
%    Can be used to index edges and faces. Given that the number of nodes 
%    in a mesh is referred to as noNum, we have the following examples
%
%    1) Edges: size_ipt = noNum*[1 1]
%    2) Triangular faces: size_ipt = noNum*[1 1 1]
%    3) Quadrilateral faces: size_ipt = noNum*[1 1 1 1]

% Define the name ['gbl___' name_ipt] as a global variable
str = sprintf('global gbl___%s', name_ipt);
eval(str)

% Store its name
str = sprintf('gbl___%s.name = ''%s'';', name_ipt, name_ipt);
eval(str)

% Store the dimension and size of the database
%   size_ipt(1) = number of elements in dimension #1
%   size_ipt(2) = number of elements in dimension #2
%   ...
%   size_ipt(n) = number of elements in dimension #n
str = sprintf('gbl___%s.size = [%s];', name_ipt, num2str(size_ipt));
eval(str)

% Store an empty data set
str = sprintf('gbl___%s.data = [];', name_ipt);
eval(str)