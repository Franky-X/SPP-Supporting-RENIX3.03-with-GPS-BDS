function [date, epoch, type, station, data, index, rcvpos, time] = readrinexobs303(r_o_name)

% 打开文件
fileID = fopen(r_o_name, 'r');
station = r_o_name(1:end-4);
disp(station)
% 检查文件是否成功打开
if fileID == -1
    error('无法打开文件');
end

epoch_cnt = -1;
epoch_buf = [];
data_buf = [];
index_buf = [];
time_buf = [];
read_state = -1;
% sys_cnt = 0;
% sysCell = {};
% typeCell = {};

% 循环读取每一行
while ~feof(fileID)
    line = fgets(fileID); % 读取一行
    if contains(line, 'TIME OF FIRST OBS')
        % 提取日期信息
        dateStr = strtrim(line(1:end-35)); % 假设日期始终位于相同的位置

        date = str2num(dateStr);
    end
    if contains(line, 'APPROX POSITION XYZ')
        rcvpos = str2num(line(1:end-22));
    end
    if contains(line, 'SYS / # / OBS TYPES')
        % 提取日期信息
        sysStr = strtrim(line(1)); % 假设日期始终位于相同的位置
        numStr = strtrim(line(6));
        num = str2double(numStr);
        typeStr = strtrim(line(7:7+num*4)); % 假设日期始终位于相同的位置
        cellArray = strsplit(typeStr, ' ');
        if sysStr == 'C'  % C for BDS, G for GPS
            type = cellArray';
        end
    end
    if contains(line, '>')
        read_state = 1;
        time_buf(end+1,:) = str2num(line(2:end));
        if epoch_cnt == -1
            epoch_cnt = epoch_cnt + 1
        else
            epoch_cnt = epoch_cnt + 30; %v2 30
        end
    end
    if read_state == 1 && ~contains(line, '>')
        sysStr = strtrim(line(1));
        if sysStr == 'C' 
            epoch_buf(end+1) = epoch_cnt;
            index_buf(end+1) = str2double(line(2:3));
   
            tmp = nan(1,9);
            tmp_data = str2num(line(4:end));
            tmp(1:length(tmp_data)) = tmp_data;
            data_buf(end+1,:) = tmp;
        end
    end
end
epoch = epoch_buf';
index = index_buf';
data = data_buf;
time = time_buf;
fclose(fileID); % 关闭文件
return;
end

% fclose(fileID); % 关闭文件
% return;