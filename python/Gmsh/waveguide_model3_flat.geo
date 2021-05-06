// Gmsh project created on Tue Mar 23 21:34:04 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {0.1, -0.05, 0, 1.0};
//+
Point(2) = {0.1, 0.05, 0, 1.0};
//+
Point(3) = {-0.1, 0.05, 0, 1.0};
//+
Point(4) = {-0.1, -0.05, 0, 1.0};
//+
Line(5) = {3, 2};
//+
Line(6) = {2, 1};
//+
Line(7) = {1, 4};
//+
Line(8) = {4, 3};
//+
Curve Loop(1) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 0.4} {
  Surface{1};
}//+
Physical Volume("freeSpace", 17) = {1};
//+
Physical Surface("bound", 18) = {2};
//+
Physical Surface("bound", 18) += {3};
//+
Physical Surface("bound", 18) += {4};
//+
Physical Surface("bound", 18) += {5};
//+
Physical Surface("port1", 19) = {1};
//+
Physical Surface("port2", 20) = {6};