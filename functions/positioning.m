function positioning(obs,nav,doy,S_path)
% Calculate user positioning by LSE solution
% Inputs: 
%        obs     = observation data
%        nav     = navigation data
%        doy     = date of year
%        S_path  = results path
%
% Save files:
%       userpos  = user position
%       disterr  = positioning error (m)
%       model    = atmospheric model (Tropo and Iono delay)
%       refpos   = reference position

%% Setting#2
stt     = 0;    % start UTC time(hour)
stp     = 25;   % stop  UTC time(hour)
int     = 30;   % interval second
% stt     = 0;    % start UTC time(hour)
% stp     = 0.25;   % stop  UTC time(hour)
% int     = 1;   % interval second
len_matrix = (stp - stt) * 3600 - 1;
% ======================================
% Ref position (Read PPP from PPPindex.txt)
refpos = reffromIPPindex(obs.station);
disp(['Calculate positioning at ' obs.station ' station'])
if isempty(refpos)
    try
        refpos = obs.rcvpos;
    catch
    end
end
refpos_lla = xyz2lla(refpos(1),refpos(2),refpos(3));
% Constants
gpscons

%% 1. Prepare matrix
% NaN
userpos.llh          = nan(len_matrix,3);  % Lat Long height
userpos.xyz          = nan(len_matrix,3);  % ECEF X Y Z
userpos.bs           = nan(len_matrix,32); % Satellite clock bias
userpos.br           = nan(len_matrix,1);  % Receiver clock bias
userpos.elevation    = nan(len_matrix,32); % elevation angle
disterr.horizontal   = nan(len_matrix,1);  % Positioning error: Horizontal
disterr.EW           = nan(len_matrix,1);  % Positioning error: East-West
disterr.NS           = nan(len_matrix,1);  % Positioning error: North-South
disterr.height       = nan(len_matrix,1);  % Positioning error: Heigth
model.tropo          = nan(len_matrix,32); % Tropospheric model
model.iono           = nan(len_matrix,32); % Ionospheric model
model.hdop           = nan(len_matrix,1); % HDOP

G = zeros(4,4);

for mode = 1:3 % correction mode (For loop#0)
% Initial USER position
x0   = zeros(4,1);
xyz0 = x0(1:3);    % initial position
br   = x0(4);      % initial bias
% mode 1 = No atmospheric delay correction
% mode 2 = Ionospheric  delay correction
% mode 3 = Tropospheric delay correction
% mode 4 = Tropospheric+Ionospheric delay correction 
    for i = ((stt*3600)+1):int:((stp*3600)-1) % receiver time (UTC)(For loop#1)

        eph_t = find(obs.epoch == (i-1));     % Check time        
        PRN   = obs.index(eph_t);             % Check satellite number (PRN)

        if length(PRN)< 4; continue;end
        
        % Second Of Day
        SOD = obs.epoch(eph_t)+((obs.date(4)*60*60)+(obs.date(5)*60)+obs.date(6));

        %% 2. Read pseudorange from Observation file
        % Pseudorange
        C1 = obs.data(eph_t,ismember(obs.type,'C2I'));

    
        %% 3. Calculate satellite position and satellite clock bias (For loop#2)
        satpos        = nan(3,length(PRN));
        satclock      = nan(length(PRN),1);
        for ii = 1:length(PRN)
            [satpos(1:3,ii),satclock(ii)] = satpos_xyz_sbias(SOD(ii),...
                PRN(ii),nav.eph,nav.index,obs.date,C1(ii));
        end
        
        if length(satclock(~isnan(satclock))) < 4;continue;end % Check number of satellite
    

        %% 4. Calculate xyz position
        % reset matrix every loop (For loop#1)
        dIon_klob     = nan(length(PRN),1);
        delay.tropo   = nan(length(PRN),1);
        delay.iono    = nan(length(PRN),1);

        % initial dx
        dx   = x0 +inf; % initial differential xyz
        stop = 0;
        
        while norm(dx) > 10^-4 % Difference of positioning
            stop = stop+1; if stop==30;break;end % iterative divergent
            % Calculate elevation angle
            [elev,~] = calelevation(satpos,xyz0);
            userlla  = xyz2lla(xyz0(1),xyz0(2),xyz0(3));
                   
            % Troposheric delay (Tropo model)
            [d_hyd,d_wet] = EstimateTropDalay(userlla(1),userlla(3),str2double(doy));
            delay.tropo = (d_hyd + d_wet).*(1.001./sqrt(0.002001 + sind(elev).^2));
          
            % Ionospheric delay (Klobuchar model) (For loop#3)
            for ii = 1:length(PRN) 
                dIon_klob(ii) = klobuchar_model(userlla(1),userlla(2),elev(ii),SOD(ii),nav.ionprm); % nano sec  newklob nav.ionprm
            end
            delay.iono = c.*dIon_klob;  % delay from Klobuchar model

            % elevation angle cutoff 15 degree
            if stop > 1
                elev_m          = elevcut(elev',elev',elev_mask);
                C1_m            = elevcut(C1',elev',elev_mask);
                satpos_m        = elevcut(satpos,elev',elev_mask);
                sclockbias_m    = elevcut(satclock',elev',elev_mask);
                d_tropo         = elevcut(delay.tropo',elev',elev_mask);
                d_iono          = elevcut(delay.iono',elev',elev_mask);
                PRN_m           = elevcut(PRN',elev,elev_mask);
            else
                delay.tropo     = 0.*elev;
                delay.iono      = 0.*elev;
                elev_m          = elevcut(elev',elev',0);
                C1_m            = elevcut(C1',elev',0);
                satpos_m        = elevcut(satpos,elev',0);
                sclockbias_m    = elevcut(satclock',elev',0);
                d_tropo         = elevcut(delay.tropo',elev',0);
                d_iono          = elevcut(delay.iono',elev',0);
                PRN_m           = elevcut(PRN',elev,0);
            end
            
            cP1 =   nan(1,length(C1_m));     
            % Pseudorange correction
            if length(cP1)~=length(sclockbias_m);break;end % Check number of satellite
            if mode == 1 % No atmospheric delay correction
                cP1 =   C1_m + c.*sclockbias_m;
            end
            if mode == 2 % Ionospheric delay correction
                cP1 =   C1_m + c.*sclockbias_m - d_iono;
            end
            if mode == 3 % Tropospheric delay correction
                cP1 =   C1_m + c.*sclockbias_m - d_tropo;
            end
            if mode == 4 % Tropospheric+Ionospheric delay correction
                cP1 =   C1_m + c.*sclockbias_m - d_tropo - d_iono;
            end
            % remove outlier
            ctest = sqrt( sum((refpos(:)*ones(1,length(satpos_m),1) - satpos_m).^2)); %Ref distance
            cError = cP1-ctest;
            cP1(abs(cError)>1000) = [];
            satpos_m(:,abs(cError)>1000) = [];
            if length(satpos_m) < 4 || length(C1_m) < 4 || length(C1_m)~=length(satpos_m);break;end % Check number of satellite
            
            length(C1_m)
            % Satellite Flight Correction with bc and clock error (SACNAG)     
            satpos_cor = nan(3,length(satpos_m));
    
            for iii = 1:length(satpos_m)
                dtflight = (cP1(iii)-br)/c + sclockbias_m(iii); % + d_iono/c
                satpos_cor(:,iii) = FlightTimeCorr(satpos_m(:,iii),dtflight); % REF
            end   
   
            % Calculate line of sight vectors and ranges from satellite to xo
            v = xyz0(:)*ones(1,length(satpos_cor),1) - satpos_cor;
            range = sqrt( sum(v.^2) );
            v = v./(ones(3,1)*range); % line of sight unit vectors from sv to xo
        
            % Calculate the a-priori range residual
            prHat = range(:) + br;
            P = cP1' - prHat;                    % matrix P
            H = [v',ones(length(satpos_cor),1)]; % matrix H
    
            % LSE solution    
            dx = inv(H'*H)*H'*P;
    
            % update xyz0 and br
            xyz0 = xyz0(:)+dx(1:3);
            br   = br+dx(4);
            G = inv(H'*H);
    
        end

        
        HDOP = sqrt(G(1,1) + G(2,2));
        
        userpos.hdop(i)            = HDOP;            % USER HDOP
        userpos.xyz(i,1:3)         = xyz0;            % USER position in xyz
        userpos.llh(i,1:3)         = xyz2lla(xyz0(1),xyz0(2),xyz0(3)); % USER position in LLA
        userpos.bs(i,PRN_m)        = sclockbias_m;    % Satellite clock bias (sec)
        userpos.br(i)              = br;              % Receiver clock bias (m)
        userpos.elevation(i,PRN_m) = elev_m;          % elevation angle (degree)
        model.tropo(i,PRN_m)       = d_tropo;         % Tropo delay (m)
        model.iono(i,PRN_m)        = d_iono;          % Iono delay  (m)
        
        if isempty(refpos)
            refpos     = nanmedian(userpos.xyz);
            refpos_lla = nanmedian(userpos.llh); 
        end
        
        % positioning error
        [disterr.horizontal(i),~]   = lldistm([userlla(1) userlla(2)],[refpos_lla(1) refpos_lla(2)]);           % Horizontal
        [disterr.EW(i),~]           = lldistm([userlla(1) refpos_lla(2)],[refpos_lla(1) refpos_lla(2)]);        % EW
        [disterr.NS(i),~]           = lldistm([refpos_lla(1) userlla(2)],[refpos_lla(1) refpos_lla(2)]);        % NS
        disterr.height(i)           = userlla(3)-refpos_lla(3);                                                 % Vertical

        numtime = 070526 + i; % 示例数值
        %base模式要除以30
        numtime = obs.time(obs.epoch(eph_t(1))/30+1,6) + obs.time(obs.epoch(eph_t(1))/30+1,5)*100 + obs.time(obs.epoch(eph_t(1))/30+1,4)*10000;
        utc_time = sprintf('%.2f', numtime); % 转换为字符串，保留两位小数
        latitude =  userpos.llh(i,1); % 31度01.4767846分
        longitude = userpos.llh(i,2); % 121度26.3320568分
        lat_dir = 'N';
        long_dir = 'E';
        pos_fix_indicator = 1; % SPP定位
        satellites_used = length(C1_m);
        hdop = HDOP;
        altitude = userpos.llh(i,3);
        altitude_units = 'M';
        geoid_separation = 11.094;
        geoid_units = 'M';
        age_of_diff_corr = '0.0';
        diff_ref_station_id = '0000';

        gngga_sentence = createGNGGA(utc_time, latitude, lat_dir, longitude, long_dir, ...
            pos_fix_indicator, satellites_used, hdop, altitude, altitude_units, ...
            geoid_separation, geoid_units, age_of_diff_corr, diff_ref_station_id);

        status = 'A'; % A=有效定位
        lat_dir = 'N';
        long_dir = 'E';
        speed = 0.00; % 地面速度
        course = 0.00; % 地面航向
        date = '261123'; % UTC日期
        variation = 0.0; % 磁偏角
        var_dir = 'E'; % 磁偏角方向
        mode_indicator = 'A'; % 模式指示符，A=自主定位

        gnrmc_sentence = createGNRMC(utc_time, status, latitude, lat_dir, longitude, long_dir, speed, course, date, variation, var_dir, mode_indicator);
%         % 打印GNGGA句式
%         disp(gngga_sentence)

        % 定义要写入的文件名
        filename = 'GNGGA_sentence_BaseStationV2_BDS.txt';

        % 打开文件准备附加写入，'a'表示附加模式
        fileID = fopen(filename, 'a');

        % 检查文件是否成功打开
        if fileID == -1
            error('File cannot be opened for writing.');
        end

        % 将GNGGA句式附加写入文件
        fprintf(fileID, '%s\r\n', gnrmc_sentence);
        fprintf(fileID, '%s\r\n', gngga_sentence);

        % 关闭文件
        fclose(fileID);

        % 显示信息表示操作完成
        disp(['GNGGA sentence appended to ', filename]);

  
    end
    
    %% Save file
    year  = num2str(obs.date(1));
    month = num2str(obs.date(2),'%.2d');
    date  = num2str(obs.date(3),'%.2d');
    m = num2str(mode);
    namepos1 = ['userpos_mode_' m];
    namepos2 = ['disterr_mode_' m];
    namepos3 = ['model_mode_' m];
    eval([namepos1 '= userpos;'])
    eval([namepos2 '= disterr;'])
    eval([namepos3 '= model;'])
    filename = [S_path 'mode_' m ...
        '\pos_m' m '_' obs.station '_' year '_' month '_' date];
    save(filename,namepos1,namepos2,namepos3,'refpos')
    disp(['Save results ... ' num2str(mode*25) '... percent'])
end
    disp([obs.station ' station Complete'])
end

    
