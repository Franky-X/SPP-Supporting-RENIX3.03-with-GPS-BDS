function gngga_sentence = createGNGGA(utc_time, latitude, lat_dir, longitude, long_dir, pos_fix_indicator, satellites_used, hdop, altitude, altitude_units, geoid_separation, geoid_units, age_of_diff_corr, diff_ref_station_id)
    % 创建GNGGA句式的函数
    % 参数列表：
    % utc_time: UTC时间，格式为'hhmmss.ss'
    % latitude: 纬度，十进制度格式
    % lat_dir: 纬度方向，'N' 或 'S'
    % longitude: 经度，十进制度格式
    % long_dir: 经度方向，'E' 或 'W'
    % pos_fix_indicator: 定位质量指标
    % satellites_used: 使用的卫星数量
    % hdop: 水平精度因子
    % altitude: 海拔高度，单位米
    % altitude_units: 海拔高度单位
    % geoid_separation: 大地水准面高度
    % geoid_units: 大地水准面高度单位
    % age_of_diff_corr: 差分GPS数据年龄
    % diff_ref_station_id: 差分参考站ID
    
    % 将十进制度纬度转换为度、分、秒的格式
    lat_deg = floor(latitude);
    lat_min = (latitude - lat_deg) * 60;
    latitude_formatted = sprintf('%02d%011.8f', lat_deg, lat_min);
    
    % 将十进制度经度转换为度、分、秒的格式
    long_deg = floor(longitude);
    long_min = (longitude - long_deg) * 60;
    longitude_formatted = sprintf('%03d%011.8f', long_deg, long_min);
    
    % 组装GNGGA句式
    gngga_sentence = sprintf('$GNGGA,%s,%s,%s,%s,%s,%d,%02d,%1.1f,%f,%s,%f,%s,%s,%s', ...
                             utc_time, latitude_formatted, lat_dir, longitude_formatted, ...
                             long_dir, pos_fix_indicator, satellites_used, hdop, ...
                             altitude, altitude_units, geoid_separation, geoid_units, ...
                             age_of_diff_corr, diff_ref_station_id);
    
    % 计算校验和
    checksum = calculateNMEAChecksum(gngga_sentence);
    
    % 添加校验和到GNGGA句式
    gngga_sentence = [gngga_sentence, '*', checksum];
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
