function marker_points = readmark(filepath)
    % READMARK : 从标记文件 (.hcm) 中读取标记点
    % 输入:
    %   filepath - 标记点文件路径
    % 输出:
    %   marker_points - [x, t] 格式的标记点坐标

    % 打开文件
    fid = fopen(filepath, 'r');
    if fid == -1
        error('无法打开文件：%s', filepath);
    end

    try
        % 读取第一行标题，确保格式符合要求
        header = fgetl(fid);
        if ~contains(header, 'Mark')
            error('文件格式错误，标题行应包含 "Mark"');
        end
        
        % 使用 textscan 读取数据部分，制表符分隔
        data = textscan(fid, '%f%f%f', 'Delimiter', '\t');
        
        % 关闭文件
        fclose(fid);

        % 提取前两列数据 (x, t)，忽略第三列
        marker_points = [data{1}, data{2}];
    catch err
        fclose(fid);
        rethrow(err); % 重新抛出错误信息
    end

    % 验证是否读取了有效数据
    if isempty(marker_points) || size(marker_points, 2) ~= 2
        error('未找到有效标记点数据！');
    end
end
