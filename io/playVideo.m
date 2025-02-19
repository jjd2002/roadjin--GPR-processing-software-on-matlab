function playVideo()
    % 获取当前脚本所在的目录（确保相对路径能够正确指向视频）
    currentDir = fileparts(mfilename('fullpath'));  % 获取当前脚本的完整路径
    videoPath = fullfile(currentDir, 'videos', '1.mp4');  % 使用相对路径指向视频

    % 检查视频文件是否存在
    if exist(videoPath, 'file') ~= 2
        errordlg('视频文件不存在，请检查路径！', '文件错误');
        return;
    end

    % 创建视频读取对象
    try
        v = VideoReader(videoPath);
    catch
        errordlg('无法读取视频文件，请检查视频格式！', '读取错误');
        return;
    end
    
    % 创建一个新的 figure 用于显示视频
    figure('Name', '视频播放', 'NumberTitle', 'off', 'MenuBar', 'none');
    ax = axes('Position', [0.1 0.3 0.8 0.6]);  % 视频显示区域
    
    % 播放视频
    while hasFrame(v)
        % 获取当前视频帧
        frame = readFrame(v);
        
        % 显示当前帧
        imshow(frame, 'Parent', ax);
        
        % 控制播放的帧率
        pause(1 / v.FrameRate);
    end
end

