cl1 = 1;
scale=4000;
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 0.3*scale, 0, 1.0};
Point(3) = {0, 0.4*scale, 0, 1.0};
Point(4) = {0, 1*scale, 0, 1.0};

Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Extrude {1*scale, 0, 0} {
  Line{1, 2, 3};
}
Transfinite Line {1, 5, 4, 6, 10, 12, 14, 3} = 8 Using Bump 1;
Transfinite Line {8, 2} = 4 Using Bump 1;

Physical Surface(50) = {7, 15};
Physical Surface(60) = {11};
Coherence;
