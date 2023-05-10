MinX=-100;
MaxX=200;
MinY=0;
MaxY=25;
R=0.5;
RefineMinX=-7;
RefineMaxX=20;
RefineMaxY=5.0;
Refine2MinX=-1.5;
Refine2MaxX=15;
Refine2MaxY=1.5;
globalRefinementFactor=2;
RefinementFactorCorarse=2;
RefinementFactorFine=0.1;
RefinementFactorFinest=0.05;
dxCoarse=RefinementFactorCorarse*globalRefinementFactor;
dxFine=RefinementFactorFine*globalRefinementFactor;
dxFinest=RefinementFactorFinest*globalRefinementFactor;
 
 
Point(1) = {MinX, MinY, 0, dxCoarse};
Point(2) = {MinX, MaxY, 0, dxCoarse};
Point(3) = {MaxX, MaxY, 0, dxCoarse};
Point(4) = {MaxX, MinY, 0, dxCoarse};
Point(5) = {-R,MinY,0,dxFinest};
Point(6) = {0,R,0,dxFinest};
Point(7) = {R,MinY,0,dxFinest};
Point(8) = {0,0,0,dxFinest};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,7};
Circle(5) = {7,8,6};
Circle(6) = {6,8,5};
Line(7) = {5,1};
 
Line Loop(21) = {7,1,2,3,4,5,6};
 
Plane Surface(22) = {21};
 
Physical Curve('inlet') = 1;
Physical Curve('symmetry') = {4,7};
Physical Curve('top') = 2;
Physical Curve('outlet') = 3;
Physical Curve('cylinder') = {5,6};
Physical Surface('interior') = 22;

