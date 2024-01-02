function [d1m,d2m]=lldistm(latlon1,latlon2)
% Calculate diferential position
% Inputs:
%   latlon1: latlon of origin point [lat lon]
%   latlon2: latlon of destination point [lat lon]
%
% Outputs:
%   d1m: distance calculated by Haversine formula
%   d2m: distance calculated based on Pythagoran theorem
%
% Example:
%   latlon1=[13.7 100];
%   latlon2=[13.7 101];
%   [d1km d2km]=distance(latlon1,latlon2)
Re = 6371.009*10^3;

lat1=latlon1(1)*pi/180;
lat2=latlon2(1)*pi/180;

lon1=latlon1(2)*pi/180;
lon2=latlon2(2)*pi/180;

deltaLat=lat2-lat1;
deltaLon=lon2-lon1;

a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2;
c=2*atan2(sqrt(a),sqrt(1-a));

d1m=Re*c;               % Haversine distance

x=deltaLon*cos((lat1+lat2)/2);
y=deltaLat;
d2m=Re*sqrt(x*x + y*y); % Pythagoran distance

end