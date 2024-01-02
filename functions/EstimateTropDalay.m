function [d_hyd,d_wet] = EstimateTropDalay(Latitude,Height,DOY)
% === Meteorologycal Parameter for Tropospheric Delay ===
% [d_hyd,d_wet] = EstimateTropDalay(Latitude,Height,DOY)
% Latitude = Latitude of a receiver on the earth
% Height = Altitude of a receiver on the earth
% DOY = day of year
% By Somkit Sophan (18 March 2019)
% -------------------------------------------------------
if isnan(Latitude)|| isnan(DOY)
    error('Value of Latitude and DOY are NaN');
end
Lat = abs(Latitude);
% == RTCA DO-229D (A.4.2.4)==
if Latitude>=0
    Dmin = 28; % for northern latitude
else
    Dmin = 211; % for southern latitude
end
if Height >=1000;Height = 1000;end % clear outage
H = Height;
k1 = 77.604; % K/mbar
k2 = 382000; % K^2/mbar
Rd = 287.054; % J/(kg.K)
gm = 9.784; % m/s^2
g = 9.80665; % m/s^2
if Lat<=15
    P0 = 1013.25;   % Pressure(mbar)
    T0 = 299.65;    % Temperature(K)
    e0 = 26.31;     % Water vapor presure(mbar)
    B0 = 6.30e-3;   % Temperature lapse rate(K/m)
    Lam0 = 2.77;    % Water vapor "lapse rate" (dimensionless)
    Delta_P = 0;
    Delta_T = 0;
    Delta_e = 0;
    Delta_B = 0;
    Delta_Lam = 0;
elseif Lat>15 && Lat<=30
    P0 = Interpol(Lat,30,15,1017.25,1013.25);
    T0 = Interpol(Lat,30,15,294.15,299.65);
    e0 = Interpol(Lat,30,15,21.79,26.31);
    B0 = Interpol(Lat,30,15,6.05e-3,6.3e-3);
    Lam0 = Interpol(Lat,30,15,3.15,2.77);
    Delta_P = Interpol(Lat,30,15,-3.75,0);
    Delta_T = Interpol(Lat,30,15,7,0);
    Delta_e = Interpol(Lat,30,15,8.85,0);
    Delta_B = Interpol(Lat,30,15,0.25e-3,0);
    Delta_Lam = Interpol(Lat,30,15,0.33,0);
elseif Lat>30 && Lat<=45
    P0 = Interpol(Lat,45,30,1015.75,1017.25);
    T0 = Interpol(Lat,45,30,283.15,294.15);
    e0 = Interpol(Lat,45,30,11.66,21.79);
    B0 = Interpol(Lat,45,30,5.58e-3,6.05e-5);
    Lam0 = Interpol(Lat,45,30,2.57,3.15);
    Delta_P = Interpol(Lat,45,30,-2.25,-3.75);
    Delta_T = Interpol(Lat,45,30,11,7);
    Delta_e = Interpol(Lat,45,30,7.24,8.85);
    Delta_B = Interpol(Lat,45,30,0.32e-3,0.25e-3);
    Delta_Lam = Interpol(Lat,45,30,0.46,0.33);
elseif Lat>45 && Lat<=60
    P0 = Interpol(Lat,60,45,1011.75,1015.75);
    T0 = Interpol(Lat,60,45,272.15,283.15);
    e0 = Interpol(Lat,60,45,6.78,11.66);
    B0 = Interpol(Lat,60,45,5.39e-3,5.58e-3);
    Lam0 = Interpol(Lat,60,45,1.81,2.57);
    Delta_P = Interpol(Lat,60,45,-1.75,-2.25);
    Delta_T = Interpol(Lat,60,45,15,11);
    Delta_e = Interpol(Lat,60,45,5.36,7.24);
    Delta_B = Interpol(Lat,60,45,0.81e-3,0.32e-3);
    Delta_Lam = Interpol(Lat,60,45,0.74,0.46);
elseif Lat>60 && Lat<=75
    P0 = Interpol(Lat,60,75,1013,1011.75);
    T0 = Interpol(Lat,60,75,263.65,272.15);
    e0 = Interpol(Lat,60,75,4.11,6.78);
    B0 = Interpol(Lat,60,75,4.53e-3,5.39e-3);
    Lam0 = Interpol(Lat,60,75,1.55,1.81);
    Delta_P = Interpol(Lat,60,75,-0.5,-1.75);
    Delta_T = Interpol(Lat,60,75,14.5,15);
    Delta_e = Interpol(Lat,60,75,3.39,5.36);
    Delta_B = Interpol(Lat,60,75,0.62e-3,0.81e-3);
    Delta_Lam = Interpol(Lat,60,75,0.30,0.74);
elseif Lat>70
    P0 = 1013;
    T0 = 263.65;
    e0 = 4.11;
    B0 = 4.53e-3;
    Lam0 = 1.55;
    Delta_P = -0.5;
    Delta_T = 14.5;
    Delta_e = 3.39;
    Delta_B = 0.62e-3;
    Delta_Lam = 0.30;
end
P = P0 - Delta_P*cos(2*pi*(DOY-Dmin)/365.25);
T = T0 - Delta_T*cos(2*pi*(DOY-Dmin)/365.25);
e = e0 - Delta_e*cos(2*pi*(DOY-Dmin)/365.25);
B = B0 - Delta_B*cos(2*pi*(DOY-Dmin)/365.25);
Lam = Lam0 - Delta_Lam*cos(2*pi*(DOY-Dmin)/365.25);
Zhyd = ((10^-6)*k1*Rd*P)/gm;
Zwet = (((10^-6)*k2*Rd)/(gm*(Lam+1)-B0*Rd))*(e/T);
d_hyd = ((1-(B*H/T))^(g/(Rd*B)))*Zhyd;
d_wet = ((1-(B*H/T))^(((Lam+1)*g/(Rd*B))-1))*Zwet;
end
