// Gmsh project created on Mon Mar 15 15:42:18 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {-0.1, 0.1, 0, 1.0};
//+
Point(2) = {0.1, 0.1, 0, 1.0};
//+
Point(3) = {0.1, -0.1, 0, 1.0};
//+
Point(4) = {-0.1, -0.1, 0, 1.0};
//+
Point(5) = {-0.09, 0.09, 0, 1.0};
//+
Point(6) = {0.09, 0.09, 0, 1.0};
//+
Point(7) = {0.09, -0.09, 0, 1.0};
//+
Point(8) = {-0.09, -0.09, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 6};
//+
Line(3) = {6, 5};
//+
Line(4) = {5, 1};
//+
Line(5) = {2, 3};
//+
Line(6) = {7, 3};
//+
Line(7) = {7, 6};
//+
Line(8) = {7, 8};
//+
Line(9) = {8, 4};
//+
Line(10) = {4, 3};
//+
Line(11) = {8, 5};
//+
Line(12) = {1, 4};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, -7, 6, -5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 9, 10, -6};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {11, 4, 12, -9};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {3, -11, -8, 7};
//+
Plane Surface(5) = {5};
//+
Extrude {0, 0, 0.4} {
  Surface{1}; 
}
//+
Extrude {0, 0, 0.4} {
  Surface{2}; 
}
//+
Extrude {0, 0, 0.4} {
  Surface{3}; 
}
//+
Extrude {0, 0, 0.4} {
  Surface{4}; 
}
//+
Extrude {0, 0, 0.4} {
  Surface{5}; 
}
//+
Physical Volume("freeSpace", 53) = {5};
//+
Physical Volume("metal", 54) = {1};
//+
Physical Volume("metal", 54) += {4};
//+
Physical Volume("metal", 54) += {3};
//+
Physical Volume("metal", 54) += {2};
//+
Physical Surface("port1", 55) = {30};
//+
Physical Surface("port2", 56) = {5};
//+
Physical Surface("test", 57) = {12, 21, 16, 8};
//+
Physical Surface("bound", 58) = {1, 2, 3, 4, 10, 15, 20, 25};
//+
Physical Surface("bound", 58) += {6};
//+
Physical Surface("bound", 58) += {14};
//+
Physical Surface("bound", 58) += {21, 18};
//+
Physical Surface(" bound", 58) -= {21};
//+
Physical Surface("bound", 58) += {23};
