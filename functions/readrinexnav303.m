function [eph, index, ionprm] = readrinexnav303(r_n_name)

% 打开文件
fileID = fopen(r_n_name, 'r');

% 检查文件是否成功打开
if fileID == -1
    error('无法打开文件');
end
ionprm = zeros(2,4);
readstate = 0;
index_buf = [];
eph_buf = [];
line_cnt = 0;
eph_cnt = 0;
% 循环读取每一行
while ~feof(fileID)
    line = fgets(fileID); % 读取一行
    line_cnt = line_cnt + 1;
    if line_cnt > 5
        if contains(line, 'C')
            readstate = 1;
            index_buf(end+1) = str2double(line(2:3));
            tmp_data = str2num(line(4:end));
        end
        if readstate == 1
           eph_cnt = eph_cnt + 1;
           tmp_data = [tmp_data str2num(line)];
           if eph_cnt == 8
               eph_buf(:,end+1) = tmp_data;
               readstate = 0;
               eph_cnt = 0;
           end
        end
    end
end
index = index_buf;
eph = eph_buf';
% disp(eph_buf)
fclose(fileID); % 关闭文件
return;
end

% fclose(fileID); % 关闭文件
% return;