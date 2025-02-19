function do_readmark(IPD)
    % DO_READMARK: 读取标记文件并在当前数据剖面上显示标记点

    % 选择标记文件
    [filename, pathname] = uigetfile('*.hcm', '选择标记文件');
    if isequal(filename, 0)
        disp('没有选择文件.');
        return;
    end
    
    % 读取标记文件内容
    filepath = fullfile(pathname, filename);
    marker_points = readmark(filepath);  % 读取标记点

    % 获取横坐标 x 和纵坐标 tt2w
    x = IPD.x;  % 横坐标
    t = IPD.tt2w;  % 纵坐标

    % 将标记点第一列除以 100，对应横坐标 x
    marker_x = marker_points(:, 1) / 100;  % 将第一列坐标转换为与 IPD.x 对应的值

    % 查找当前已存在的图形窗口
    datafig = findobj('tag', 'datafigure');
    procdatafig = findobj('tag', 'procdatafigure');  % 查找处理后的数据窗口

    % 如果没有找到窗口，显示错误
    if isempty(datafig) && isempty(procdatafig)
        errordlg('没有找到数据图形窗口。');
        return;
    end

    % 根据当前显示的图形窗口选择绘制标记的目标窗口
    if ~isempty(datafig)
        % 在 datafig 窗口中绘制标记
        figure(datafig);
        hold on;
        plot(marker_x, zeros(size(marker_x)), 'r+', 'MarkerSize', 10);  % 绘制红色十字形标记点
        hold off;
        title('添加标记点 - 输入数据');
    elseif ~isempty(procdatafig)
        % 在 procdatafig 窗口中绘制标记
        figure(procdatafig);
        hold on;
        plot(marker_x, zeros(size(marker_x)), 'r+', 'MarkerSize', 10);  % 绘制红色十字形标记点
        hold off;
        title('添加标记点 - 处理数据');
    end
end

function marker_points = readmark(filepath)
    % READMARK: 读取标记文件
    % 读取文件并返回标记点（以矩阵形式）

    fid = fopen(filepath, 'r');
    if fid == -1
        error('无法打开文件: %s', filepath);
    end

    % 跳过文件开头的标头部分
    fgetl(fid);  % 跳过第一行 "Mark"

    % 读取标记点
    marker_points = [];
    while ~feof(fid)
        line = fgetl(fid);  % 逐行读取
        if ischar(line)
            % 将每行数据分隔并转为数值
            tokens = str2double(strsplit(line));
            if numel(tokens) >= 3  % 只处理有3列的行
                marker_points = [marker_points; tokens(1:3)];  % 存储前3列
            end
        end
    end

    fclose(fid);
end

