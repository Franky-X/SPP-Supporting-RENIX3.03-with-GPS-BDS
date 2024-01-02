function [Out] = Interpol(Wp,Wmax,Wmin,Vmax,Vmin)
% ======== Interpolation ========================
% [Out] = Interpol(Wp,Wmax,Wmin,Vmax,Vmin)
% Wp = weighting point
% Wmax = Max weight
% Wmin = Min weight
% Vmax = Value of max weight
% Vmin = Value of min weight
% By Somkit Sophan (18 March 2019)
% ===== RTCA DO-229D (A.4.2.4)=====
Out = Vmin + (Vmax-Vmin)*((Wp-Wmin)/(Wmax-Wmin));
end