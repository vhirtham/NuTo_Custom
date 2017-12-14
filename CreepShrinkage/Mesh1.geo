W = 0.2;
H = 0.2;
R = 0.015;
R1 = 0.02;
R2 = 0.03;
X1 = W * 0.65;
Y1 = H * 0.7;
X2 = W * 0.45;
Y2 = H * 0.4;


SB = (W+H) / 2 /40;
SI = (W+H) / 2 /100;

Point(1)  = {   0.0,     0.0,   0.0,   SB};
Point(2)  = {   0.0,       H,   0.0,   SB};
Point(3)  = {     W,       H,   0.0,   SB};
Point(4)  = {     W,     0.0,   0.0,   SB};
Point(5)  = {    X1,      Y1,   0.0,   SI};
Point(6)  = {    X1,   Y1+R1,   0.0,   SI};
Point(7)  = {    X1,   Y1-R1,   0.0,   SI};
Point(8)  = {    X2,      Y2,   0.0,   SI};
Point(9)  = {    X2,   Y2+R2,   0.0,   SI};
Point(10) = {    X2,   Y2-R2,   0.0,   SI};


Line(1)    = {1,2};
Line(2)    = {2,3};
Line(3)    = {3,4};
Line(4)    = {4,1};

Circle(5)  = {6, 5, 7};
Circle(6)  = {7, 5, 6};
Circle(7)  = {9, 8, 10};
Circle(8)  = {10, 8, 9};

Line Loop(15) = {1,2,3,4};
Line Loop(16) = {5,6};
Line Loop(17) = {7,8};

Plane Surface(21) = {15,16,17};

Plane Surface(22) = {16};
Plane Surface(23) = {17};

Physical Surface(27) = {21};
Physical Surface(28) = {22,23};
