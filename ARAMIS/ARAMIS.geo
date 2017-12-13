W = 0.2;
H = 0.2;
R = 0.015;
L = 0.047;
CX = W/2;
CY = H/2;

SB = (W+H) / 2 /40;
SI = (W+H) / 2 /100;

Point(1)  = {   0.0,     0.0,   0.0,   SB};
Point(2)  = {   0.0,       H,   0.0,   SB};
Point(3)  = {     W,       H,   0.0,   SB};
Point(4)  = {     W,     0.0,   0.0,   SB};
Point(5)  = {    CX,      CY,   0.0,   SI};
Point(6)  = {    CX,    CY+R,   0.0,   SI};
Point(7)  = {    CX,    CY-R,   0.0,   SI};
Point(8)  = {  CX+L,    CY-L,   0.0,   SI};
Point(9)  = {  CX+L,  CY-L+R,   0.0,   SI};
Point(10) = {  CX+L,  CY-L-R,   0.0,   SI};
Point(11) = {  CX-L,    CY+L,   0.0,   SI};
Point(12) = {  CX-L,  CY+L+R,   0.0,   SI};
Point(13) = {  CX-L,  CY+L-R,   0.0,   SI};
Point(14) = {  CX+L,    CY+L,   0.0,   SI};
Point(15) = {  CX+L,  CY+L+R,   0.0,   SI};
Point(16) = {  CX+L,  CY+L-R,   0.0,   SI};
Point(17) = {  CX-L,    CY-L,   0.0,   SI};
Point(18) = {  CX-L,  CY-L+R,   0.0,   SI};
Point(19) = {  CX-L,  CY-L-R,   0.0,   SI};


Line(1)    = {1,2};
Line(2)    = {2,3};
Line(3)    = {3,4};
Line(4)    = {4,1};

Circle(5)  = {6, 5, 7};
Circle(6)  = {7, 5, 6};
Circle(7)  = {9, 8, 10};
Circle(8)  = {10, 8, 9};
Circle(9)  = {12, 11, 13};
Circle(10) = {13, 11, 12};
Circle(11) = {15, 14, 16};
Circle(12) = {16, 14, 15};
Circle(13) = {18, 17, 19};
Circle(14) = {19, 17, 18};

Line Loop(15) = {1,2,3,4};
Line Loop(16) = {5,6};
Line Loop(17) = {7,8};
Line Loop(18) = {9,10};
Line Loop(19) = {11,12};
Line Loop(20) = {13,14};

Plane Surface(21) = {15,16,17,18,19,20};

Plane Surface(22) = {16};
Plane Surface(23) = {17};
Plane Surface(24) = {18};
Plane Surface(25) = {19};
Plane Surface(26) = {20};


Physical Surface(27) = {21};
Physical Surface(28) = {24, 25, 22, 23, 26};
