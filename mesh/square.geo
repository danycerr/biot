cl1 = 1;
scale=4000;
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 1, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {1, 0, 0, 1.0};
Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Dilate {{scale, scale, scale}, 1} {
  Surface{6};
}
Transfinite Line {1, 4, 3, 2} = 4 Using Bump 1;
Dilate {{0, 0, 0}, 4000} {
  Surface{6};
}
