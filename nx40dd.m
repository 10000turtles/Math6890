%Filewrittenbyheat1dImp.C
xa=0;xb=1;kappa=0.1;t=1;maxErr= 4.229e-16;cpuTimeStep= 4.210e-04;
Nx=40;dx=  2.500000e-02;numGhost=1;n1a=0;n1b=40;nd1a=-1;nd1b=41;
solutionName='polyDD';
x=[-2.500000000000000e-020.000000000000000e+002.500000000000000e-025.000000000000000e-027.500000000000001e-021.000000000000000e-011.250000000000000e-011.500000000000000e-01...
1.750000000000000e-012.000000000000000e-012.250000000000000e-012.500000000000000e-012.750000000000000e-013.000000000000000e-013.250000000000000e-013.500000000000000e-01...
3.750000000000000e-014.000000000000000e-014.250000000000000e-014.500000000000000e-014.750000000000000e-015.000000000000000e-015.250000000000000e-015.500000000000000e-01...
5.750000000000001e-016.000000000000001e-016.250000000000000e-016.500000000000000e-016.750000000000000e-017.000000000000001e-017.250000000000001e-017.500000000000000e-01...
7.750000000000000e-018.000000000000000e-018.250000000000001e-018.500000000000001e-018.750000000000000e-019.000000000000000e-019.250000000000000e-019.500000000000001e-01...
9.750000000000001e-011.000000000000000e+001.025000000000000e+00];
u=[1.185187499999999e+001.200000000000000e+001.215187500000000e+001.230750000000000e+001.246687500000000e+001.263000000000000e+001.279687500000000e+001.296750000000000e+00...
1.314187500000000e+001.332000000000000e+001.350187500000000e+001.368750000000000e+001.387687500000000e+001.407000000000000e+001.426687500000000e+001.446750000000000e+00...
1.467187499999999e+001.487999999999999e+001.509187500000000e+001.530749999999999e+001.552687499999999e+001.574999999999999e+001.597687499999999e+001.620749999999999e+00...
1.644187499999999e+001.667999999999999e+001.692187499999999e+001.716749999999999e+001.741687500000000e+001.766999999999999e+001.792687500000000e+001.818750000000000e+00...
1.845187500000000e+001.872000000000000e+001.899187500000000e+001.926750000000000e+001.954687500000000e+001.983000000000000e+002.011687500000000e+002.040750000000000e+00...
2.070187500000000e+002.100000000000000e+002.130187500000000e+00];
err=[0.000000000000000e+000.000000000000000e+000.000000000000000e+000.000000000000000e+000.000000000000000e+000.000000000000000e+002.220446049250313e-164.440892098500626e-16...
0.000000000000000e+00-2.220446049250313e-160.000000000000000e+000.000000000000000e+00-2.220446049250313e-16-4.440892098500626e-16-4.440892098500626e-16-2.220446049250313e-16...
-4.440892098500626e-16-6.661338147750939e-16-4.440892098500626e-16-6.661338147750939e-16-6.661338147750939e-16-6.661338147750939e-16-8.881784197001252e-16-6.661338147750939e-16...
-6.661338147750939e-16-8.881784197001252e-16-6.661338147750939e-16-6.661338147750939e-16-4.440892098500626e-16-6.661338147750939e-16-6.661338147750939e-16-2.220446049250313e-16...
0.000000000000000e+000.000000000000000e+00-2.220446049250313e-16-2.220446049250313e-160.000000000000000e+002.220446049250313e-160.000000000000000e+000.000000000000000e+00...
4.440892098500626e-160.000000000000000e+001.175156250000000e-01];
figure(1)
plot(x,u)
title("Solution")
x_label("x")
y_label("u")
figure(2)
plot(x,err)
title("Error")
x_label("x")
y_label("err")
