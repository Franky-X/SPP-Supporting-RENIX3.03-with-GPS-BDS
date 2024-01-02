function y=f(xi,yi,zi,vx,vy,vz,xls,yls,zls)
% format long e;
GM=398600.44E9;%-m3/s2
ae = 6378136; %地球半径；单位m
J0_2= 1082625.7E-9;%地球重力系数；
we = 7.292115e-5;%地球自转速度；
r = sqrt(xi*xi+yi*yi+zi*zi);
f1 = -GM*xi/r^3-1.5*J0_2*GM*ae*ae*xi*(1-5*zi^2/r^2)/r^5+we*we*xi+2*we*vy+xls;
f2 = -GM*yi/r^3-1.5*J0_2*GM*ae*ae*yi*(1-5*zi^2/r^2)/r^5+we*we*yi-2*we*vx+yls;
f3 = -GM*zi/r^3-1.5*J0_2*GM*ae*ae*zi*(3-5*zi^2/r^2)/r^5+zls;
y = [f1,f2,f3];
end

