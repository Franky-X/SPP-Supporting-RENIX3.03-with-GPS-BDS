function newsatpos = FlightTimeCorr(satpos, dTflightSeconds)
% To correct satellite position from pseudorange
% Inputs:
%         satpos          = Satellite position
%         dTflightSeconds = Time error (second)
% Output:
%         newsatpos       = Corrected satellite position

we = 7.2921151467e-5;           %   Earth rotation rate (rad/sec)

theta = we * dTflightSeconds;

R3 = [ cos(theta)    sin(theta)   0;
      -sin(theta)    cos(theta)   0;
       0                0         1];
   
newsatpos = (R3*satpos(:))';