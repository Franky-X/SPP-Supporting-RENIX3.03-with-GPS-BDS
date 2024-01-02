r_o_name = 'BaseStation.obs'; % observation file's name
r_n_name = 'BaseStation.nav'; % navigation file's name (if blank, program will downloaded from IGS)
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
[obs_test, nav_test, doy_test,Year_test] = readrinex303(r_o_name,r_n_name,R_path); 
