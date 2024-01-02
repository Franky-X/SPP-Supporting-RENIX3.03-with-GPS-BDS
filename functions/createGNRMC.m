function gnrmc_sentence = createGNRMC(utc_time, status, latitude, lat_dir, longitude, long_dir, speed, course, date, variation, var_dir, mode_indicator)
    % 创建GNRMC句式的函数
    % 参数列表：
    % utc_time: UTC时间，格式为'hhmmss.ss'
    % status: 状态，A=有效定位，V=无效定位
    % latitude: 纬度，十进制度格式
    % lat_dir: 纬度方向，'N' 或 'S'
    % longitude: 经度，十进制度格式
    % long_dir: 经度方向，'E' 或 'W'
    % speed: 地面速度，单位节
    % course: 地面航向，度
    % date: UTC日期，格式为'ddmmyy'
    % variation: 磁偏角
    % var_dir: 磁偏角方向，'E' 或 'W'
    % mode_indicator: 模式指示符（A=自主定位，D=差分，E=估计，N=数据无效）
    
    % 将十进制度纬度转换为度和分的格式
    lat_deg = floor(latitude);
    lat_min = (latitude - lat_deg) * 60;
    latitude_formatted = sprintf('%02d%010.7f', lat_deg, lat_min);
    
    % 将十进制度经度转换为度和分的格式
    long_deg = floor(longitude);
    long_min = (longitude - long_deg) * 60;
    longitude_formatted = sprintf('%03d%010.7f', long_deg, long_min);
    
    % 组装GNRMC句式
    gnrmc_sentence = sprintf('$GNRMC,%s,%s,%s,%s,%s,%s,%.2f,%.2f,%s,%.1f,%s,%s', ...
                             utc_time, status, latitude_formatted, lat_dir, longitude_formatted, ...
                             long_dir, speed, course, date, variation, var_dir, mode_indicator);
    
    % 计算校验和
    checksum = calculateNMEAChecksum(gnrmc_sentence);
    
    % 添加校验和到GNRMC句式
    gnrmc_sentence = [gnrmc_sentence, '*', checksum];
end

function checksum = calculateNMEAChecksum(sentence)
    % 计算NMEA句式的校验和
    checksum = 0;
    for i = 2:length(sentence) % 从第一个字符'$'后面开始
        checksum = bitxor(checksum, double(sentence(i)));
    end
    checksum = dec2hex(checksum, 2);
    if length(checksum) < 2
        checksum = ['0', checksum]; % 确保校验和为两位数
    end
end
