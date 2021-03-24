// Gmsh project created on Tue Mar 23 21:34:04 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {0.1, -0.1, 0, 1.0};
//+
Point(2) = {0.1, 0.1, 0, 1.0};
//+
Point(3) = {-0.1, 0.1, 0, 1.0};
//+
Point(4) = {-0.1, -0.1, 0, 1.0};
//+
Point(5) = {-0.09, -0.09, 0, 1.0};
//+
Point(6) = {0.09, 0.09, 0, 1.0};
//+
Point(7) = {-0.09, 0.09, 0, 1.0};
//+
Point(8) = {0.09, -0.09, 0, 1.0};
//+
Line(1) = {7, 6};
//+
Line(2) = {6, 8};
//+
Line(3) = {8, 5};
//+
Line(4) = {7, 5};
//+
Line(5) = {3, 2};
//+
Line(6) = {2, 1};
//+
Line(7) = {1, 4};
//+
Line(8) = {4, 3};
//+
Curve Loop(1) = {8, 5, 6, 7};
//+
Curve Loop(2) = {1, 2, 3, -4};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {1, 2, 3, -4};
//+
Plane Surface(2) = {3};
//+
Extrude {0, 0, 0.4} {
  Surface{2}; Surface{1}; 
}
//+
Physical Volume("freeSpace", 53) = {2};
//+
Physical Volume("metal", 54) = {1};//+
Physical Surface("port1", 55) = {7};
//+
Physical Surface("port2", 56) = {2};
//+
Physical Surface("bound", 57) = {9, 8, 10, 11, 12, 1};
