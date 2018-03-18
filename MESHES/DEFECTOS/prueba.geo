//parameters
sc = 5; //separacion central
sl = 3; //separacion lateral
h = 2; //altura defectos
b = 2; //base defectos
H = 8*h; //altura
B = sc+2*sl+2*b; //base
y = 1*H/3; //posicion "y" de defectos

lado_ele = 0.05;

// points
//contorno
Point(1) = {0, 0, 0, lado_ele};
Point(2) = {B, 0, 0, lado_ele};
Point(3) = {B, H, 0, lado_ele};
Point(4) = {0, H, 0, lado_ele};
//defecto izquierdo
Point(5) = {sl, y-h/2, 0, lado_ele};
Point(6) = {sl+b, y-h/2, 0, lado_ele};
Point(7) = {sl+b, y+h/2, 0, lado_ele};
Point(8) = {sl, y+h/2, 0, lado_ele};
//defecto derecho
Point(9) = {B-sl-b, y-h/2, 0, lado_ele};
Point(10) = {B-sl, y-h/2, 0, lado_ele};
Point(11) = {B-sl, y+h/2, 0, lado_ele};
Point(12) = {B-sl-b, y+h/2, 0, lado_ele};
//linea de carga
Point(13) = {B, H/6, 0, lado_ele};
Point(14) = {0, H/6, 0, lado_ele};

//lines
//contornos
Line(1) = {1, 2};
Line(2) = {2, 13};
Line(3) = {13, 3};
Line(4) = {3, 4};
Line(5) = {4, 14};
Line(6) = {14, 1};
//defecto izquierdo
Line(7) = {5, 6};
Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 5};
//defecto derecho
Line(11) = {9, 10};
Line(12) = {10, 11};
Line(13) = {11, 12};
Line(14) = {12, 9};
//linea de carga
Line(15) = {13, 14};

//Surfaces
Line Loop(1) = {1, 2, 15, 6};
Line Loop(2) = {-15, 3, 4, 5};
Line Loop(3) = {7, 8, 9, 10};
Line Loop(4) = {11, 12, 13, 14};
Plane Surface(1) = {1}; //superficie inferior
Plane Surface(2) = {2, 3, 4}; //superficie superior

//Physical groups
//Physical Line(1) = {2, 3, 5, 6}; //laterales
Physical Line(100) = {1}; //base y tope
Physical Line(300) = {15}; //carga
Physical Surface(1000) = {1};
Physical Surface(2000) = {2};
Physical Line(500) = {4}; //respuesta