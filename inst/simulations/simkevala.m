addpath('kevala2017')
[z1hat,z2hat, eps,del] = FactorScoresEst('save/d.csv',[3,3]);
zz = dlmread('save/du.csv');
[factores,splinedeg,Knopt,Knots,z1hat0] = BSplineEst(z1hat,z2hat,1);
B = bspline_basismatrix(splinedeg+1,Knots,zz(1,:)');
yhat = B*factores;
results = [zz(1,:)', yhat, z1hat{1}, z2hat]
dlmwrite('save/dkevala.csv', results')
