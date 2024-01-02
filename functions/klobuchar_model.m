function dIon = klobuchar_model(fi,lambda,elev,tow,ionprm)
% Ionospheric delay calculation (Klobuchar model,1987)
% Inputs: 
%       fi     = latitude of rcv position
%       lambda = longitude of rcv position
%       elev   = elevation angle
%       tow    = time of week
%       ionprm = navigation coeficients (alpha and beta)
% Outputs:
%         dIon = Ionospheric delay

alfa  = ionprm(1,:); % alpha
beta  = ionprm(2,:); % Beta

azimuth  = 360;
c        =  2.99792458e8;             % speed of light
deg2semi =  1./180.;                  % degees to semisircles
semi2rad =  pi;                       % semisircles to radians
deg2rad  =  pi/180.;                  % degrees to radians

a = azimuth*deg2rad;                  % asimuth in radians
e = elev*deg2semi;                    % elevation angle in
                                      % semicircles

psi = 0.0137 / (e+0.11) - 0.022;      % Earth Centered angle

lat_i = fi*deg2semi + psi*cos(a);     % Subionospheric lat
if (lat_i > 0.416)
    lat_i = 0.416;
  elseif(lat_i < -0.416)
      lat_i = -0.416;
end

                                      % Subionospheric long
long_i = lambda*deg2semi + (psi*sin(a)/cos(lat_i*semi2rad));

                                      % Geomagnetic latitude
lat_m = lat_i + 0.064*cos((long_i-1.617)*semi2rad);

t = 4.32e4*long_i + tow;
t = mod(t,86400.);                    % Seconds of day
if (t > 86400.)
    t = t - 86400.;
end
if (t < 0.)
    t = t + 86400.;
end

sF = 1. + 16.*(0.53-e)^3;             % Slant factor

                                      % Period of model
PER = beta(1) + beta(2)*lat_m + beta(3)*lat_m^2 +beta(4)*lat_m^3;

if (PER < 72000.)
    PER = 72000.;
end

x = 2.*pi*(t-50400.)/PER ;            % Phase of the model
                                      % (Max at 14.00 =
                                      % 50400 sec local time)

                                      % Amplitud of the model
AMP = alfa(1) + alfa(2)*lat_m + alfa(3)*lat_m^2 +alfa(4)*lat_m^3;
if(AMP < 0.)
    AMP = 0.;
end

                                      % Ionospheric corr.
if(abs(x) > 1.57)
    dIon1 = sF * (5.e-9);
else
    dIon1 = sF * (5.e-9 + AMP*(1. - x*x/2. + x*x*x*x/24.));
end

dIon = dIon1; % * 10.^9





