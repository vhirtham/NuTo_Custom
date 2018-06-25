RO = 50;
RI = RO/2;
H = 300;

NEO = RO / 2;
NEI = RI / 1;

// Points ----------------------------------------------------------------

Point(1)   = {   0.,  0.,  0., NEO};
Point(2)   = {   RO,  0.,  0., NEO};
Point(3)   = {   0.,  RO,  0., NEO};
Point(4)   = {  -RO,  0.,  0., NEO};
Point(5)   = {   0., -RO,  0., NEO};
Point(6)   = {   RI,  0.,  0., NEI};
Point(7)   = {   0.,  RI,  0., NEI};
Point(8)   = {  -RI,  0.,  0., NEI};
Point(9)   = {   0., -RI,  0., NEI};

Point(11)  = {   0.,  0.,   H, NEO};
Point(12)  = {   RO,  0.,   H, NEO};
Point(13)  = {   0.,  RO,   H, NEO};
Point(14)  = {  -RO,  0.,   H, NEO};
Point(15)  = {   0., -RO,   H, NEO};
Point(16)  = {   RI,  0.,   H, NEI};
Point(17)  = {   0.,  RI,   H, NEI};
Point(18)  = {  -RI,  0.,   H, NEI};
Point(19)  = {   0., -RI,   H, NEI};



// Lines -----------------------------------------------------------------

Circle(1)   = {2, 1, 3};
Circle(2)   = {3, 1, 4};
Circle(3)   = {4, 1, 5};
Circle(4)   = {5, 1, 2};
Line(5)     = {6,7};
Line(6)     = {7,8};
Line(7)     = {8,9};
Line(8)     = {9,6};
Line(9)     = {6,2};
Line(10)    = {7,3};
Line(11)    = {8,4};
Line(12)    = {9,5};

Circle(21)  = {12, 11, 13};
Circle(22)  = {13, 11, 14};
Circle(23)  = {14, 11, 15};
Circle(24)  = {15, 11, 12};
Line(25)    = {16,17};
Line(26)    = {17,18};
Line(27)    = {18,19};
Line(28)    = {19,16};
Line(29)    = {16,12};
Line(30)    = {17,13};
Line(31)    = {18,14};
Line(32)    = {19,15};

Line(41)    = {2,12};
Line(42)    = {3,13};
Line(43)    = {4,14};
Line(44)    = {5,15};
Line(45)    = {6,16};
Line(46)    = {7,17};
Line(47)    = {8,18};
Line(48)    = {9,19};



// Surfaces --------------------------------------------------------------


// Upper Face -----------------------------------

// Outer Part ----------

Line Loop(101) = {1, -10, -5, 9};
Plane Surface(201) = {101};

Line Loop(102) = {2, -11, -6, 10};
Plane Surface(202) = {102};

Line Loop(103) = {3, -12, -7, 11};
Plane Surface(203) = {103};

Line Loop(104) = {4, -9, -8, 12};
Plane Surface(204) = {104};

// Inner Part ----------

Line Loop(105) = {5, 6, 7, 8};
Plane Surface(205) = {105};


// Lower Face -----------------------------------

// Outer Part ----------

Line Loop(111) = {21, -30, -25, 29};
Plane Surface(211) = {111};

Line Loop(112) = {22, -31, -26, 30};
Plane Surface(212) = {112};

Line Loop(113) = {23, -32, -27, 31};
Plane Surface(213) = {113};

Line Loop(114) = {24, -29, -28, 32};
Plane Surface(214) = {114};

// Inner Part ----------

Line Loop(115) = {25, 26, 27, 28};
Plane Surface(215) = {115};


// Side -----------------------------------------

Line Loop(121) = {1, 42, -21, -41};
Ruled Surface(221) = {121};

Line Loop(122) = {2, 43, -22, -42};
Ruled Surface(222) = {122};

Line Loop(123) = {3, 44, -23, -43};
Ruled Surface(223) = {123};

Line Loop(124) = {4, 41, -24, -44};
Ruled Surface(224) = {124};


// Inner Planes ---------------------------------

Line Loop(131) = {5, 46, -25, -45};
Plane Surface(231) = {131};

Line Loop(132) = {6, 47, -26, -46};
Plane Surface(232) = {132};

Line Loop(133) = {7, 48, -27, -47};
Plane Surface(233) = {133};

Line Loop(134) = {8, 45, -28, -48};
Plane Surface(234) = {134};


// Connection Planes ----------------------------

Line Loop(141) = {9, 41, -29, -45};
Plane Surface(241) = {141};

Line Loop(142) = {10, 42, -30, -46};
Plane Surface(242) = {142};

Line Loop(143) = {11, 43, -31, -47};
Plane Surface(243) = {143};

Line Loop(144) = {12, 44, -32, -48};
Plane Surface(244) = {144};



// Volumes ---------------------------------------------------------------


// Outer parts ----------------------------------

Surface Loop(301) = {201, 211, 221, 231, 241, 242};
Volume (401) = {301};

Surface Loop(302) = {202, 212, 222, 232, 242, 243};
Volume (402) = {302};

Surface Loop(303) = {203, 213, 223, 233, 243, 244};
Volume (403) = {303};

Surface Loop(304) = {204, 214, 224, 234, 244, 241};
Volume (404) = {304};


// Center ---------------------------------------

Surface Loop(305) = {231, 232, 233, 234, 205, 215};
Volume (405) = {305};



// Physical Groups -------------------------------------------------------

Physical Surface("UpperFace") = {201, 202, 203, 204, 205};
Physical Surface("LowerFace") = {211, 212, 213, 214, 215};
Physical Surface("SideFace") = {221, 222, 223, 224};
Physical Volume("Volume") = {401, 402, 403, 404, 405};

