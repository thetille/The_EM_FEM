// Gmsh project created on Tue May 18 20:29:57 2021
SetFactory("OpenCASCADE");
//+
a = DefineNumber[ 0.0106684, Name "Parameters/a" ];
//+
b = DefineNumber[ 0.004318, Name "Parameters/b" ];
//+
l = DefineNumber[ 0.020, Name "Parameters/l" ];
//+
dev_space = DefineNumber[ 0.005, Name "Parameters/dev_space" ];
//+
middle = DefineNumber[ 0.014, Name "Parameters/middle" ];
//+
tria = DefineNumber[ 0.004, Name "Parameters/tria" ];
//+
Point(1) = {a/2, 0, -l, 1.0};
//+
Point(2) = {-a/2, 0, -l, 1.0};
//+
Point(3) = {-a/2, 0, 0, 1.0};
//+
Point(4) = {-a-(dev_space/2), 0, middle, 1.0};
//+
Point(5) = {-a-(dev_space/2), 0, middle+l, 1.0};
//+
Point(6) = {-(dev_space/2), 0, middle+l, 1.0};
//+
Point(7) = {-dev_space/2, 0, middle, 1.0};
//+
Point(8) = {0, 0, middle-tria, 1.0};
//+
Point(9) = {dev_space/2, 0, middle, 1.0};
//+
Point(10) = {dev_space/2, 0, middle+l, 1.0};
//+
Point(11) = {(dev_space/2)+a, 0, middle+l, 1.0};
//+
Point(12) = {(dev_space/2)+a, 0, middle, 1.0};
//+
Point(13) = {a/2, 0, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 10};
//+
Line(9) = {10, 11};
//+
Line(10) = {11, 12};
//+
Line(11) = {12, 13};
//+
Line(12) = {13, 1};
//+
Line(13) = {1, 2};
//+
Curve Loop(1) = {12, 13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
//+
Plane Surface(1) = {1};
//+
Extrude {0, b, 0} {
  Surface{1}; 
}
//+
Physical Volume("freeSpace", 40) = {1};
//+
Physical Surface("bound", 41) = {1, 4, 5, 6, 13, 14, 2, 15, 11, 10, 9, 8};
//+
Physical Surface("port1", 42) = {3};
//+
Physical Surface("port2", 43) = {7};
//+
Physical Surface("port3", 44) = {12};
