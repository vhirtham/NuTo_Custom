W = 200;
H = 200;
R = 15;
L = 47;
CX = W/2;
CY = H/2;

SB = (W+H) / 2 /60;
SI = (W+H) / 2 /60;
SIC = (W+H) / 2 /2;

Point(1)  = {   0.0,     0.0,   0.0,   SB};
Point(2)  = {   0.0,       H,   0.0,   SB};
Point(3)  = {     W,       H,   0.0,   SB};
Point(4)  = {     W,     0.0,   0.0,   SB};
Point(5)  = {    CX,      CY,   0.0,   SIC};
Point(6)  = {    CX,    CY+R,   0.0,   SI};
Point(7)  = {    CX,    CY-R,   0.0,   SI};
Point(8)  = {  CX+L,    CY-L,   0.0,   SIC};
Point(9)  = {  CX+L,  CY-L+R,   0.0,   SI};
Point(10) = {  CX+L,  CY-L-R,   0.0,   SI};
Point(11) = {  CX-L,    CY+L,   0.0,   SIC};
Point(12) = {  CX-L,  CY+L+R,   0.0,   SI};
Point(13) = {  CX-L,  CY+L-R,   0.0,   SI};
Point(14) = {  CX+L,    CY+L,   0.0,   SIC};
Point(15) = {  CX+L,  CY+L+R,   0.0,   SI};
Point(16) = {  CX+L,  CY+L-R,   0.0,   SI};
Point(17) = {  CX-L,    CY-L,   0.0,   SIC};
Point(18) = {  CX-L,  CY-L+R,   0.0,   SI};
Point(19) = {  CX-L,  CY-L-R,   0.0,   SI};




Line(1)    = {1,2};
Line(2)    = {2,3};
Line(3)    = {3,4};
Line(4)    = {4,1};
Line(13)    = {5,6};
Line(14)    = {5,7};
Line(15)    = {8,9};
Line(16)    = {8,10};
Line(17)    = {11,12};
Line(18)    = {11,13};
Line(19)    = {14,15};
Line(20)    = {14,16};
Line(21)    = {17,18};
Line(22)    = {17,19};

Circle(43)  = {6, 5, 7};
Circle(44)  = {7, 5, 6};
Circle(45)  = {9, 8, 10};
Circle(46)  = {10, 8, 9};
Circle(47)  = {12, 11, 13};
Circle(48) = {13, 11, 12};
Circle(49) = {15, 14, 16};
Circle(50) = {16, 14, 15};
Circle(51) = {18, 17, 19};
Circle(52) = {19, 17, 18};


// Edges -----------------------------------------------------------------

// Outer Edges

Line Loop(63) = {1,2,3,4};


// Central cylinder

Line Loop(69) = {13,43,-14};
Line Loop(70) = {14,44,-13};


// Right lower cylinder
Line Loop(76) = {15,45,-16};
Line Loop(77) = {16,46,-15};

// Left upper cylinder
Line Loop(83) = {17,47,-18};
Line Loop(84) = {18,48,-17};

// Right upper cylinder
Line Loop(90) = {19,49,-20};
Line Loop(91) = {20,50,-19};


// Left lower cylinder
Line Loop(97) = {21,51,-22};
Line Loop(98) = {22,52,-21};

// Helper
Line Loop(104) = {43,44};
Line Loop(106) = {45,46};
Line Loop(108) = {47,48};
Line Loop(110) = {49,50};
Line Loop(112) = {51,52};



// Surfaces --------------------------------------------------------------

// Central cylinder
Plane Surface(204) = {69};
Plane Surface(205) = {70};

// Right lower cylinder
Plane Surface(211) = {76};
Plane Surface(212) = {77};

// Left upper cylinder
Plane Surface(218) = {83};
Plane Surface(219) = {84};

// Right upper cylinder
Plane Surface(225) = {90};
Plane Surface(226) = {91};

// Left lower cylinder
Plane Surface(232) = {97};
Plane Surface(233) = {98};

// Front and Back
Plane Surface(239) = {63,104,106,108,110,112};





// Physical Groups -------------------------------------------------------

//Physical Surface("LeftMatrixBoundary") = {200};
//Physical Surface("RightMatrixBoundary") = {202};
//Physical Surface("TopMatrixBoundary") = {201};
//Physical Surface("BottomMatrixBoundary") = {203};
//Physical Surface("FrontMatrixBoundary") = {239};
//Physical Surface("BackMatrixBoundary") = {240};

Physical Surface("Matrix") = {-239};
Physical Surface("Granite") = {204,205,211,212,218,219,225,226,232,233};

