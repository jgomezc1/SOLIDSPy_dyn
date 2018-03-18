// Parameters
L = 6.0;
H1 = 3.0;
H2 = 3.0;
hc = L/10.0;
// Points
Point(1) = {0.0, 0.0, 0, hc};
Point(2) = {L, 0.0, 0, hc};
Point(3) = {L, H1, 0, hc};
Point(4) = {L, H1+H2, 0, hc};
Point(5) = {L/2, H1+H2, 0, hc};
Point(6) = {0, H1+H2, 0, hc};
Point(7) = {0, H1, 0, hc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};
Line(8) = {3, 7};

// Surfaces
Line Loop(8) = {1, 2, 8, 7};
Plane Surface(9) = {8};
Line Loop(10) = {-8, 3, 4, 5 , 6};
Plane Surface(11) = {10};

// Physical groups
Physical Surface(100) = {9};  // Top layer
Physical Surface(200) = {11};  // Bottom layr
Physical Line(400) = {1};  // Bottom line
Physical Point(500) = {5};  // Load point

