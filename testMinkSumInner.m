clear
clc
A = Polyhedron('lb',-ones(4,1),'ub',ones(4,1));
B = A;
C = A + B;
D = Polyhedron('V',allcomb([-1,0,1],[-1,0,1],[-1,0,1],[-1,0,1]));
disp('>>> Dinner = minkSumInner(A, B, D)');
tic
Dinner = minkSumInner(A, B, D);
toc;
C.contains(Dinner);

disp('>>> Einner = minkSumInner(A, B)');
tic
Einner = minkSumInner(A, B);
toc;
C.contains(Einner);

