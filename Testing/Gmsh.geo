Mesh.MshFileVersion = 2.2;lc=1.0;//Characteristic length

Point(0)={-40.1906200200441, 107.189905102499, 0, lc};
Point(1)={67.1755776941667, 107.189905102499, 0, lc};
Point(2)={-40.1906200200441, 17.8159745449888, 0, lc};
Point(3)={67.1755776941667, 17.8159745449888, 0, lc};

Point(4)={-44.1906200200441, 107.189905102499, 0, lc};	//$NAME L1_left
Point(5)={-44.1906200200441, 17.8159745449888, 0, lc};	//$NAME L1_right
Point(6)={71.1755776941667, 17.8159745449888, 0, lc};	//$NAME L3_left
Point(7)={71.1755776941667, 107.189905102499, 0, lc};	//$NAME L3_right

Line(1)={1,0};
Line(2)={0,2};
Line(3)={2,3};
Line(4)={3,1};
Line(5)={2,5};
Line(6)={5,4};
Line(7)={4,0};
Line(8)={1,7};
Line(9)={7,6};
Line(10)={6,3};

Transfinite Line (1)=108.0 Using Progression 1.0;
Transfinite Line (2)=90.0 Using Progression 1.0;
Transfinite Line (3)=108.0 Using Progression 1.0;
Transfinite Line (4)=90.0 Using Progression 1.0;
Transfinite Line (5)=5.0 Using Progression 1.0;
Transfinite Line (6)=90.0 Using Progression 1.0;
Transfinite Line (7)=5.0 Using Progression 1.0;
Transfinite Line (8)=5.0 Using Progression 1.0;
Transfinite Line (9)=90.0 Using Progression 1.0;
Transfinite Line (10)=5.0 Using Progression 1.0;

Line Loop(1)={1,2,3,4};
Line Loop(2)={2, 5, 6, 7};
Line Loop(3)={4, 8, 9, 10};

Plane Surface(1)={1};
Plane Surface(2)={2};
Plane Surface(3)={3};

Transfinite Surface(1)={0,1,2,3};
Transfinite Surface(2)={0,2,5,4};
Transfinite Surface(3)={3,1,7,6};

Recombine Surface{1};
Recombine Surface{2};
Recombine Surface{3};
