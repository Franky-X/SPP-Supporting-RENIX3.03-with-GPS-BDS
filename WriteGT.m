p_path = [pwd '\'];             % Program path
R_path = [p_path 'RINEX\'];     % RINEX path
S_path = [p_path 'Results\'];    % Results path
path(path,[p_path 'function']);


stt     = 0;    % start UTC time(hour)
stp     = 25;   % stop  UTC time(hour)
int     = 30;   % interval second
for i = ((stt*3600)+1):int:((stp*3600)-1)
    numtime = 070526+i; % 示例数值
    utc_time = sprintf('%.2f', numtime); % 转换为字符串，保留两位小数
    latitude =  31.025026951; % 31度01.4767846分
    longitude = 121.439449785; % 121度26.3320568分
    lat_dir = 'N';
    long_dir = 'E';
    pos_fix_indicator = 4; % 假设为差分固定定位
    satellites_used = 19;
    hdop = 1.0;
    altitude = 37.0669;
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
    filename = 'BaseStationGT.txt';

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
