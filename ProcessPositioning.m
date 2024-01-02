% Positioning calculation from RINEX 3.03
% Calculate user position by using the single point method
% Original by Chaoran Xiong
% Version 1.00 
% (12/01/2024) - Create the program
% Output
% mode 1 = No atmospheric delay correction
% mode 2 = Ionospheric  delay correction
% mode 3 = Tropospheric delay correction
% mode 4 = Tropospheric+Ionospheric delay correction

close all;clear
warning off
tic

% RINEX file
r_o_name = 'BaseStation-v2.obs'; % observation file's name
r_n_name = 'BaseStation.nav'; % navigation file's name (if blank, program will downloaded from IGS)
% r_o_name = 'Rover.obs'; % observation file's name
% r_n_name = 'Rover.nav'; % navigation file's name (if blank, program will downloaded from IGS)
% r_n_name = ''; % navigation file's name

% Setting#1
% =========== Program's path ==========================
p_path = [pwd '\'];             % Program path
R_path = [p_path 'RINEX\'];     % RINEX path
S_path = [p_path 'Results\'];    % Results path
path(path,[p_path 'function']);

% 1. Read RINEX (using readrinex .mex file)
% Check file
checkfileRN(r_o_name,R_path);
% Read RINEX  

%  [obs,nav,doy,Year] = readrinex211(r_o_name,r_n_name,R_path); 
[obs,nav,doy,Year] = readrinex303(r_o_name,r_n_name,R_path); 
year  = num2str(obs.date(1));
month = num2str(obs.date(2),'%.2d');
date  = num2str(obs.date(3),'%.2d');

% 2. Estimate user positioning
positioning(obs,nav,doy,S_path);


% 3. Plot Error
ploterror(year,month,date,obs.station,S_path)
toc
% remove file (reset)
% for m = 1:4; rm = num2str(m);delete([S_path 'mode_' rm '\*_' obs.station '_' year '_' month '_' date '.mat']);end
warning on

