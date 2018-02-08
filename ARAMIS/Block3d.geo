W = 200;
H = 200;
D = 100;
R = 15;
L = 47;
CX = W/2;
CY = H/2;

SB = (W+H) / 2 /20;
SI = (W+H) / 2 /40;
SIC = (W+H) / 2 /10;

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
Point(21) = {   0.0,     0.0,     D,   SB};
Point(22) = {   0.0,       H,     D,   SB};
Point(23) = {     W,       H,     D,   SB};
Point(24) = {     W,     0.0,     D,   SB};
Point(25) = {    CX,      CY,     D,   SIC};
Point(26) = {    CX,    CY+R,     D,   SI};
Point(27) = {    CX,    CY-R,     D,   SI};
Point(28) = {  CX+L,    CY-L,     D,   SIC};
Point(29) = {  CX+L,  CY-L+R,     D,   SI};
Point(30) = {  CX+L,  CY-L-R,     D,   SI};
Point(31) = {  CX-L,    CY+L,     D,   SIC};
Point(32) = {  CX-L,  CY+L+R,     D,   SI};
Point(33) = {  CX-L,  CY+L-R,     D,   SI};
Point(34) = {  CX+L,    CY+L,     D,   SIC};
Point(35) = {  CX+L,  CY+L+R,     D,   SI};
Point(36) = {  CX+L,  CY+L-R,     D,   SI};
Point(37) = {  CX-L,    CY-L,     D,   SIC};
Point(38) = {  CX-L,  CY-L+R,     D,   SI};
Point(39) = {  CX-L,  CY-L-R,     D,   SI};



Line(1)    = {1,2};
Line(2)    = {2,3};
Line(3)    = {3,4};
Line(4)    = {4,1};
Line(5)    = {1,21};
Line(6)    = {2,22};
Line(7)    = {3,23};
Line(8)    = {4,24};
Line(9)    = {21,22};
Line(10)    = {22,23};
Line(11)    = {23,24};
Line(12)    = {24,21};
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
Line(23)    = {25,26};
Line(24)    = {25,27};
Line(25)    = {28,29};
Line(26)    = {28,30};
Line(27)    = {31,32};
Line(28)    = {31,33};
Line(29)    = {34,35};
Line(30)    = {34,36};
Line(31)    = {37,38};
Line(32)    = {37,39};
Line(33)    = {6,26};
Line(34)    = {7,27};
Line(35)    = {9,29};
Line(36)    = {10,30};
Line(37)    = {12,32};
Line(38)    = {13,33};
Line(39)    = {15,35};
Line(40)    = {16,36};
Line(41)    = {18,38};
Line(42)    = {19,39};

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
Circle(53)  = {26, 25, 27};
Circle(54)  = {27, 25, 26};
Circle(55)  = {29, 28, 30};
Circle(56)  = {30, 28, 29};
Circle(57)  = {32, 31, 33};
Circle(58) = {33, 31, 32};
Circle(59) = {35, 34, 36};
Circle(60) = {36, 34, 35};
Circle(61) = {38, 37, 39};
Circle(62) = {39, 37, 38};


// Edges -----------------------------------------------------------------

// Outer Edges

Line Loop(63) = {1,2,3,4};
Line Loop(64) = {9,10,11,12};
Line Loop(65) = {1,6,-9,-5};
Line Loop(66) = {2,7,-10,-6};
Line Loop(67) = {3,8,-11,-7};
Line Loop(68) = {4,5,-12,-8};


// Central cylinder

Line Loop(69) = {13,43,-14};
Line Loop(70) = {14,44,-13};
Line Loop(71) = {23,53,-24};
Line Loop(72) = {24,54,-23};
Line Loop(73) = {13,33,-23,24,-34,-14};
Line Loop(74) = {43,34,-53,-33};
Line Loop(75) = {44,33,-54,-34};


// Right lower cylinder
Line Loop(76) = {15,45,-16};
Line Loop(77) = {16,46,-15};
Line Loop(78) = {25,55,-26};
Line Loop(79) = {26,56,-25};
Line Loop(80) = {15,35,-25,26,-36,-16};
Line Loop(81) = {45,36,-55,-35};
Line Loop(82) = {46,35,-56,-36};

// Left upper cylinder
Line Loop(83) = {17,47,-18};
Line Loop(84) = {18,48,-17};
Line Loop(85) = {27,57,-28};
Line Loop(86) = {28,58,-27};
Line Loop(87) = {17,37,-27,28,-38,-18};
Line Loop(88) = {47,38,-57,-37};
Line Loop(89) = {48,37,-58,-38};

// Right upper cylinder
Line Loop(90) = {19,49,-20};
Line Loop(91) = {20,50,-19};
Line Loop(92) = {29,59,-30};
Line Loop(93) = {30,60,-29};
Line Loop(94) = {19,39,-29,30,-40,-20};
Line Loop(95) = {49,40,-59,-39};
Line Loop(96) = {50,39,-60,-40};


// Left lower cylinder
Line Loop(97) = {21,51,-22};
Line Loop(98) = {22,52,-21};
Line Loop(99) = {31,61,-32};
Line Loop(100) = {32,62,-31};
Line Loop(101) = {21,41,-31,32,-42,-22};
Line Loop(102) = {51,42,-61,-41};
Line Loop(103) = {52,41,-62,-42};

// Helper
Line Loop(104) = {43,44};
Line Loop(105) = {53,54};
Line Loop(106) = {45,46};
Line Loop(107) = {55,56};
Line Loop(108) = {47,48};
Line Loop(109) = {57,58};
Line Loop(110) = {49,50};
Line Loop(111) = {59,60};
Line Loop(112) = {51,52};
Line Loop(113) = {61,62};


// Surfaces --------------------------------------------------------------

// Sides
Plane Surface(200) = {65};
Plane Surface(201) = {66};
Plane Surface(202) = {67};
Plane Surface(203) = {68};




// Front and Back
Plane Surface(239) = {63};
Plane Surface(240) = {64};


// Volumes ---------------------------------------------------------------

Surface Loop(300) = {200,201,202,203,239,240};
Volume(301) = {300};


Physical Volume(312) = {301};

