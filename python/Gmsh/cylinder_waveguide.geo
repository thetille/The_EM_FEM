// Gmsh project created on Tue Mar 23 10:42:36 2021
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 0.4, 0.125, 2*Pi};
//+
Physical Volume("aire", 4) = {1};
//+
Physical Surface("port1", 5) = {2};
//+
Physical Surface("port2", 6) = {3};
//+
Physical Surface("bound", 7) = {1};
