



% function OPD = do_ttoz(IPD, v1d)
% %
% % Callback function to drive the routine "ttoz.m"
% % Time to depth conversion for constant or layered velocity structures
% %
% % Copyright (C) 2005, Andreas Tzanis. All rights reserved.
% %
% 
% OPD = discardprocdata; 
% % Trap common errors ...
% if isempty(IPD.d), 
%    erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
%    uiwait(erh);  
% end
% 
% if isempty(v1d); 
%     erh = warndlg('Please provide an 1-D velocity model!', ...
%         'MATGPR : WARNING!'); 
%     uiwait(erh); 
% end 
% 
% % Proceed with time-to-depth conversion
% OPD = IPD;
% [OPD.z, OPD.ns, OPD.dz, OPD.d] = ttoz(IPD.tt2w, IPD.d, IPD.dt, v1d, 0, 0);   
% if isempty(OPD.d)
%     disp('TTOZ > Operation aborted - No O/P data returned!');
% end
% OPD.zlab = 'Depth (m)';
% viewdata(OPD.x, OPD.z, OPD.d, 'outdata', OPD.xlab, OPD.zlab); 
% 
% iss = size(OPD.history, 1); 
% OPD.history(iss + 1, 1) = cellstr('Performed Time to Depth Conversion');
% 
% 
% 
% 
% 
% 
% % Initialize empty data table and counter for point numbering
% columnNames = {'序号', '距离起始点距离 (m)', '雷达与桩基直线距离 (m)', '实际埋深 (m)'};
% dataTable = uitable('Parent', gcf, 'Units', 'normalized', 'Position', [0.86 0.4 0.12 0.3], ...
%     'ColumnName', columnNames, 'Data', {}, 'CellEditCallback', @(src, evt)updateTable(src, evt));  % Add callback for table update
% 
% % Initialize pointCounter in the app data for sequential numbering
% setappdata(gcf, 'pointCounter', 0);  % Store pointCounter
% 
% % Add middle-click functionality for depth calculation
% set(gcf, 'WindowButtonDownFcn', @(~,~)onMouseClick(gcf, OPD.x, OPD.z, dataTable));
% 
% % Add button to export table data to Excel
% uicontrol('Style', 'pushbutton', 'String', '导出为 Excel', 'Units', 'normalized', ...
%     'Position', [0.86 0.1 0.12 0.1], 'Callback', @(src, evt)exportToExcel(dataTable));  % Add Export button
% 
% end
% 
% function onMouseClick(fig, x, z, dataTable)
%     % Get the selected point
%     clickType = get(fig, 'SelectionType');
%     if strcmp(clickType, 'extend')  % Middle-click detected
%         point = get(gca, 'CurrentPoint');
%         xClicked = point(1, 1);
%         zClicked = point(1, 2);
% 
%         % Find the closest x and z from OPD.xlab and OPD.zlab
%         [~, idxX] = min(abs(x - xClicked));  % Find closest x
%         [~, idxZ] = min(abs(z - zClicked));  % Find closest z
%         xClicked = x(idxX);  % Update x to the closest value
%         zClicked = z(idxZ);  % Update z to the closest value
% 
%         % Ask for offset input
%         prompt = {'请输入地面间隔 (m):'};
%         dlgTitle = '输入计算深度的参数';
%         numLines = 1;
%         defaultAns = {'0.5'};
%         answer = inputdlg(prompt, dlgTitle, numLines, defaultAns);
% 
%         if ~isempty(answer)
%             offset = str2double(answer{1});
% 
%             % Validate the input
%             if isnan(offset) || offset < 0
%                 errordlg('请输入有效的数字。', '输入错误');
%                 return;
%             end
% 
%             % Calculate actual depth based on z
%             actualDepth = sqrt(zClicked^2 - offset^2);  % Modify depth calculation
% 
%             % Increment point counter
%             pointCounter = getappdata(gcf, 'pointCounter') + 1;
%             setappdata(gcf, 'pointCounter', pointCounter);  % Update the pointCounter
% 
%             % Add new row of data
%             newRow = {pointCounter, xClicked, zClicked, actualDepth};  % Create new row data
%             dataTable.Data(end + 1, :) = newRow;  % Add new data to table
% 
%             % Update the plot with the new point
%             hold on;  % Ensure new points are added to the current plot
%             plot(xClicked, zClicked, 'ro', 'MarkerSize', 8);  % Plot red dot
%             hold off;
%         end
%     end
% end
% 
% % Function to manually export table to Excel via button
% function exportToExcel(dataTable)
%     % 获取表格数据
%     data = dataTable.Data;
% 
%     % 检查是否有数据可以导出
%     if isempty(data)
%         errordlg('没有数据可以导出！', '导出错误');
%         return;
%     end
% 
%     % 获取列名
%     columnNames = dataTable.ColumnName;
% 
%     % 将列名转换为列向量，以便与数据拼接
%     columnNames = columnNames(:);  % 强制转换列名为列向量
% 
%     % 确保列数与数据列一致
%     [numRows, numCols] = size(data);
%     numColsColumnNames = numel(columnNames);  % 列名的数量
% 
%     % 如果列数不一致，填充或截断数据
%     if numCols < numColsColumnNames
%         % 用空单元格填充数据列
%         for i = 1:numRows
%             data{i, numCols+1:numColsColumnNames} = {''};
%         end
%     elseif numCols > numColsColumnNames
%         % 如果数据列数大于列名，截断多余的列
%         data = data(:, 1:numColsColumnNames);
%     end
% 
%     % 合并列名与数据
%     outputData = [columnNames'; data];  % 列名和数据合并为一个cell数组
% 
%     % 设置导出文件路径
%     filename = 'data_points.xlsx';  % Excel 文件名
%     try
%         % 将列名和数据一起写入 Excel 文件
%         writecell(outputData, filename, 'Sheet', 1, 'Range', 'A1');  % 将数据写入 Excel 文件
%         msgbox('数据已成功导出为 Excel 文件！', '导出成功');  % 弹出成功提示框
%     catch ME
%         % 如果发生错误，显示详细的错误信息
%         errordlg(['导出 Excel 文件时出错：' ME.message], '导出错误');  % 弹出详细错误信息
%     end
% end
% 






% function OPD = do_ttoz(IPD, v1d)
%     % Callback function to drive the routine "ttoz.m"
%     % Time to depth conversion for constant or layered velocity structures
%     %
%     % Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%     %
% 
%     OPD = discardprocdata; 
%     % Trap common errors ...
%     if isempty(IPD.d), 
%         erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
%         uiwait(erh);  
%     end
% 
%     if isempty(v1d); 
%         erh = warndlg('Please provide an 1-D velocity model!', ...
%             'MATGPR : WARNING!'); 
%         uiwait(erh); 
%     end 
% 
%     % Proceed with time-to-depth conversion
%     OPD = IPD;
%     [OPD.z, OPD.ns, OPD.dz, OPD.d] = ttoz(IPD.tt2w, IPD.d, IPD.dt, v1d, 0, 0);   
%     if isempty(OPD.d)
%         disp('TTOZ > Operation aborted - No O/P data returned!');
%     end
%     OPD.zlab = 'Depth (m)';
%     viewdata(OPD.x, OPD.z, OPD.d, 'outdata', OPD.xlab, OPD.zlab); 
% 
%     iss = size(OPD.history, 1); 
%     OPD.history(iss + 1, 1) = cellstr('Performed Time to Depth Conversion');
% 
%     % Initialize empty data table and counter for point numbering
%     columnNames = {'序号', 'K()', '距离起始点距离 (m)', '雷达与桩基直线距离 (m)', '实际埋深 (m)'};
%     dataTable = uitable('Parent', gcf, 'Units', 'normalized', 'Position', [0.86 0.4 0.12 0.3], ...
%         'ColumnName', columnNames, 'Data', {}, 'CellEditCallback', @(src, evt)updateTable(src, evt));  % Add callback for table update
% 
%     % Initialize pointCounter and K() value
%     setappdata(gcf, 'pointCounter', 0);  % Store pointCounter
%     setappdata(gcf, 'fixedK', 147);  % Store the fixed K value (e.g., 147)
%     setappdata(gcf, 'increment', 100);  % Store increment value for K() generation
% 
%     % Add middle-click functionality for depth calculation
%     set(gcf, 'WindowButtonDownFcn', @(~,~)onMouseClick(gcf, OPD.x, OPD.z, dataTable));
% 
%     % Add button to export table data to Excel
%     uicontrol('Style', 'pushbutton', 'String', '导出为 Excel', 'Units', 'normalized', ...
%         'Position', [0.86 0.1 0.12 0.1], 'Callback', @(src, evt)exportToExcel(dataTable));  % Add Export button
% 
% end
% 
% function onMouseClick(fig, x, z, dataTable)
%     % Get the selected point
%     clickType = get(fig, 'SelectionType');
%     if strcmp(clickType, 'extend')  % Middle-click detected
%         point = get(gca, 'CurrentPoint');
%         xClicked = point(1, 1);
%         zClicked = point(1, 2);
% 
%         % Find the closest x and z from OPD.xlab and OPD.zlab
%         [~, idxX] = min(abs(x - xClicked));  % Find closest x
%         [~, idxZ] = min(abs(z - zClicked));  % Find closest z
%         xClicked = x(idxX);  % Update x to the closest value
%         zClicked = z(idxZ);  % Update z to the closest value
% 
%         % Ask for offset input
%         prompt = {'请输入地面间隔 (m):'};
%         dlgTitle = '输入计算深度的参数';
%         numLines = 1;
%         defaultAns = {'0.5'};
%         answer = inputdlg(prompt, dlgTitle, numLines, defaultAns);
% 
%         if ~isempty(answer)
%             offset = str2double(answer{1});
% 
%             % Validate the input
%             if isnan(offset) || offset < 0
%                 errordlg('请输入有效的数字。', '输入错误');
%                 return;
%             end
% 
%             % Calculate actual depth based on z
%             actualDepth = sqrt(zClicked^2 - offset^2);  % Modify depth calculation
% 
%             % Retrieve the point counter and K() values from app data
%             pointCounter = getappdata(gcf, 'pointCounter');
%             fixedK = getappdata(gcf, 'fixedK');  % The fixed K value
%             increment = getappdata(gcf, 'increment');  % The increment value for K()
% 
%             % Increment the point counter
%             pointCounter = pointCounter + 1;
%             setappdata(gcf, 'pointCounter', pointCounter);  % Update the counter
% 
%             % Generate the K() value for the new point
%             K_value = sprintf('K%d+%d', fixedK, pointCounter * increment);
% 
%             % Add new row to the table
%             newRow = {pointCounter, K_value, xClicked, zClicked, actualDepth};
%             dataTable.Data(end + 1, :) = newRow;  % Add new data to table
% 
%             % Update the plot with the new point
%             hold on;  % Ensure new points are added to the current plot
%             plot(xClicked, zClicked, 'ro', 'MarkerSize', 8);  % Plot red dot
%             hold off;
%         end
%     end
% end
% 
% % Function to manually export table to Excel via button
% function exportToExcel(dataTable)
%     % 获取表格数据
%     data = dataTable.Data;
% 
%     % 检查是否有数据可以导出
%     if isempty(data)
%         errordlg('没有数据可以导出！', '导出错误');
%         return;
%     end
% 
%     % 获取列名
%     columnNames = dataTable.ColumnName;
% 
%     % 将列名转换为列向量，以便与数据拼接
%     columnNames = columnNames(:);  % 强制转换列名为列向量
% 
%     % 确保列数与数据列一致
%     [numRows, numCols] = size(data);
%     numColsColumnNames = numel(columnNames);  % 列名的数量
% 
%     % 如果列数不一致，填充或截断数据
%     if numCols < numColsColumnNames
%         % 用空单元格填充数据列
%         for i = 1:numRows
%             data{i, numCols+1:numColsColumnNames} = {''};
%         end
%     elseif numCols > numColsColumnNames
%         % 如果数据列数大于列名，截断多余的列
%         data = data(:, 1:numColsColumnNames);
%     end
% 
%     % 合并列名与数据
%     outputData = [columnNames'; data];  % 列名和数据合并为一个cell数组
% 
%     % 设置导出文件路径
%     filename = 'data_points.xlsx';  % Excel 文件名
%     try
%         % 将列名和数据一起写入 Excel 文件
%         writecell(outputData, filename, 'Sheet', 1, 'Range', 'A1');  % 将数据写入 Excel 文件
%         msgbox('数据已成功导出为 Excel 文件！', '导出成功');  % 弹出成功提示框
%     catch ME
%         % 如果发生错误，显示详细的错误信息
%         errordlg(['导出 Excel 文件时出错：' ME.message], '导出错误');  % 弹出详细错误信息
%     end
% end








% function OPD = do_ttoz(IPD, v1d)
%     % Callback function to drive the routine "ttoz.m"
%     % Time to depth conversion for constant or layered velocity structures
%     %
%     % Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%     %
% 
%     OPD = discardprocdata; 
%     % Trap common errors ...
%     if isempty(IPD.d), 
%         erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
%         uiwait(erh);  
%     end
% 
%     if isempty(v1d); 
%         erh = warndlg('Please provide an 1-D velocity model!', ...
%             'MATGPR : WARNING!'); 
%         uiwait(erh); 
%     end 
% 
%     % Proceed with time-to-depth conversion
%     OPD = IPD;
%     [OPD.z, OPD.ns, OPD.dz, OPD.d] = ttoz(IPD.tt2w, IPD.d, IPD.dt, v1d, 0, 0);   
%     if isempty(OPD.d)
%         disp('TTOZ > Operation aborted - No O/P data returned!');
%     end
%     OPD.zlab = 'Depth (m)';
%     viewdata(OPD.x, OPD.z, OPD.d, 'outdata', OPD.xlab, OPD.zlab); 
% 
%     iss = size(OPD.history, 1); 
%     OPD.history(iss + 1, 1) = cellstr('Performed Time to Depth Conversion');
% 
%     % Initialize empty data table and counter for point numbering
%     columnNames = {'K()', '距离起始点距离 (m)', '雷达与桩基直线距离 (m)', '实际埋深 (m)'};
%     dataTable = uitable('Parent', gcf, 'Units', 'normalized', 'Position', [0.86 0.4 0.12 0.3], ...
%         'ColumnName', columnNames, 'Data', {}, 'CellEditCallback', @(src, evt)updateTable(src, evt));  % Add callback for table update
% 
%     % Initialize pointCounter and K() value
%     setappdata(gcf, 'pointCounter', 0);  % Store pointCounter
%     setappdata(gcf, 'fixedK', 147);  % Store the fixed K value (e.g., 147)
%     setappdata(gcf, 'increment', 100);  % Store increment value for K() generation
% 
%     % Add middle-click functionality for depth calculation
%     set(gcf, 'WindowButtonDownFcn', @(~,~)onMouseClick(gcf, OPD.x, OPD.z, dataTable));
% 
%     % Add button to export table data to Excel
%     uicontrol('Style', 'pushbutton', 'String', '导出为 Excel', 'Units', 'normalized', ...
%         'Position', [0.86 0.1 0.12 0.1], 'Callback', @(src, evt)exportToExcel(dataTable));  % Add Export button
% 
% end
% 
% function onMouseClick(fig, x, z, dataTable)
%     % Get the selected point
%     clickType = get(fig, 'SelectionType');
%     if strcmp(clickType, 'extend')  % Middle-click detected
%         point = get(gca, 'CurrentPoint');
%         xClicked = point(1, 1);
%         zClicked = point(1, 2);
% 
%         % Find the closest x and z from OPD.xlab and OPD.zlab
%         [~, idxX] = min(abs(x - xClicked));  % Find closest x
%         [~, idxZ] = min(abs(z - zClicked));  % Find closest z
%         xClicked = x(idxX);  % Update x to the closest value
%         zClicked = z(idxZ);  % Update z to the closest value
% 
%         % Ask for offset input
%         prompt = {'请输入地面间隔 (m):'};
%         dlgTitle = '输入计算深度的参数';
%         numLines = 1;
%         defaultAns = {'0.5'};
%         answer = inputdlg(prompt, dlgTitle, numLines, defaultAns);
% 
%         if ~isempty(answer)
%             offset = str2double(answer{1});
% 
%             % Validate the input
%             if isnan(offset) || offset < 0
%                 errordlg('请输入有效的数字。', '输入错误');
%                 return;
%             end
% 
%             % Calculate actual depth based on z
%             actualDepth = sqrt(zClicked^2 - offset^2);  % Modify depth calculation
% 
%             % Retrieve the point counter and K() values from app data
%             pointCounter = getappdata(gcf, 'pointCounter');
%             fixedK = getappdata(gcf, 'fixedK');  % The fixed K value
%             increment = getappdata(gcf, 'increment');  % The increment value for K()
% 
%             % Increment the point counter
%             pointCounter = pointCounter + 1;
%             setappdata(gcf, 'pointCounter', pointCounter);  % Update the counter
% 
%             % Generate the K() value for the new point
%             K_value = sprintf('K%d+%d', fixedK, pointCounter * increment);
% 
%             % Add new row to the table
%             newRow = {K_value, xClicked, zClicked, actualDepth};
%             dataTable.Data(end + 1, :) = newRow;  % Add new data to table
% 
%             % Update the plot with the new point
%             hold on;  % Ensure new points are added to the current plot
%             plot(xClicked, zClicked, 'ro', 'MarkerSize', 8);  % Plot red dot
%             hold off;
%         end
%     end
% end
% 
% % Function to manually export table to Excel via button
% function exportToExcel(dataTable)
%     % 获取表格数据
%     data = dataTable.Data;
% 
%     % 检查是否有数据可以导出
%     if isempty(data)
%         errordlg('没有数据可以导出！', '导出错误');
%         return;
%     end
% 
%     % 获取列名
%     columnNames = dataTable.ColumnName;
% 
%     % 将列名转换为列向量，以便与数据拼接
%     columnNames = columnNames(:);  % 强制转换列名为列向量
% 
%     % 确保列数与数据列一致
%     [numRows, numCols] = size(data);
%     numColsColumnNames = numel(columnNames);  % 列名的数量
% 
%     % 如果列数不一致，填充或截断数据
%     if numCols < numColsColumnNames
%         % 用空单元格填充数据列
%         for i = 1:numRows
%             data{i, numCols+1:numColsColumnNames} = {''};
%         end
%     elseif numCols > numColsColumnNames
%         % 如果数据列数大于列名，截断多余的列
%         data = data(:, 1:numColsColumnNames);
%     end
% 
%     % 合并列名与数据
%     outputData = [columnNames'; data];  % 列名和数据合并为一个cell数组
% 
%     % 设置导出文件路径
%     filename = 'data_points.xlsx';  % Excel 文件名
%     try
%         % 将列名和数据一起写入 Excel 文件
%         writecell(outputData, filename, 'Sheet', 1, 'Range', 'A1');  % 将数据写入 Excel 文件
%         msgbox('数据已成功导出为 Excel 文件！', '导出成功');  % 弹出成功提示框
%     catch ME
%         % 如果发生错误，显示详细的错误信息
%         errordlg(['导出 Excel 文件时出错：' ME.message], '导出错误');  % 弹出详细错误信息
%     end
% end






function OPD = do_ttoz(IPD, v1d)
%
% Callback function to drive the routine "ttoz.m"
% Time to depth conversion for constant or layered velocity structures
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors ...
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
end

if isempty(v1d); 
    erh = warndlg('Please provide an 1-D velocity model!', ...
        'MATGPR : WARNING!'); 
    uiwait(erh); 
end 

% Proceed with time-to-depth conversion
OPD = IPD;
[OPD.z, OPD.ns, OPD.dz, OPD.d] = ttoz(IPD.tt2w, IPD.d, IPD.dt, v1d, 0, 0);   
if isempty(OPD.d)
    disp('TTOZ > Operation aborted - No O/P data returned!');
end
OPD.zlab = 'Depth (m)';
viewdata(OPD.x, OPD.z, OPD.d, 'outdata', OPD.xlab, OPD.zlab); 

iss = size(OPD.history, 1); 
OPD.history(iss + 1, 1) = cellstr('Performed Time to Depth Conversion');

% Initialize empty data table and counter for point numbering
columnNames = {'序号', '距离起始点距离 (m)', '雷达与桩基直线距离 (m)', '实际埋深 (m)'};
dataTable = uitable('Parent', gcf, 'Units', 'normalized', 'Position', [0.9 0.4 0.1 0.3], ...
    'ColumnName', columnNames, 'Data', {}, 'CellEditCallback', @(src, evt)updateTable(src, evt));  % Add callback for table update

% Initialize pointCounter in the app data for sequential numbering
setappdata(gcf, 'pointCounter', 0);  % Store pointCounter

% Add middle-click functionality for depth calculation
set(gcf, 'WindowButtonDownFcn', @(~,~)onMouseClick(gcf, OPD.x, OPD.z, dataTable));

% Add button to export table data to Excel
uicontrol('Style', 'pushbutton', 'String', '导出为 Excel', 'Units', 'normalized', ...
    'Position', [0.86 0.1 0.12 0.1], 'Callback', @(src, evt)exportToExcel(dataTable));  % Add Export button

% Add button for K() input
uicontrol('Style', 'pushbutton', 'String', '输入 K()', 'Units', 'normalized', ...
    'Position', [0.86 0.2 0.12 0.1], 'Callback', @(src, evt)inputKValue(dataTable));  % Add K() Input button

end

function onMouseClick(fig, x, z, dataTable)
    % Get the selected point
    clickType = get(fig, 'SelectionType');
    if strcmp(clickType, 'extend')  % Middle-click detected
        point = get(gca, 'CurrentPoint');
        xClicked = point(1, 1);
        zClicked = point(1, 2);

        % Find the closest x and z from OPD.xlab and OPD.zlab
        [~, idxX] = min(abs(x - xClicked));  % Find closest x
        [~, idxZ] = min(abs(z - zClicked));  % Find closest z
        xClicked = x(idxX);  % Update x to the closest value
        zClicked = z(idxZ);  % Update z to the closest value

        % Ask for offset input
        prompt = {'请输入地面间隔 (m):'};
        dlgTitle = '输入计算深度的参数';
        numLines = 1;
        defaultAns = {'0.5'};
        answer = inputdlg(prompt, dlgTitle, numLines, defaultAns);

        if ~isempty(answer)
            offset = str2double(answer{1});

            % Validate the input
            if isnan(offset) || offset < 0
                errordlg('请输入有效的数字。', '输入错误');
                return;
            end

            % Calculate actual depth based on z
            actualDepth = sqrt(zClicked^2 - offset^2);  % Modify depth calculation

            % Increment point counter
            pointCounter = getappdata(gcf, 'pointCounter') + 1;
            setappdata(gcf, 'pointCounter', pointCounter);  % Update the pointCounter

            % Add new row of data
            newRow = {pointCounter, xClicked, zClicked, actualDepth};  % Create new row data
            dataTable.Data(end + 1, :) = newRow;  % Add new data to table

            % Update the plot with the new point
            hold on;  % Ensure new points are added to the current plot
            plot(xClicked, zClicked, 'ro', 'MarkerSize', 8);  % Plot red dot
            hold off;
        end
    end
end

% Function to update the table and export data to Excel
function updateTable(src, evt)
    % Get the updated data from the table
    data = src.Data;

    % You can optionally print or process data here
    disp('Table data updated:');
    disp(data);

    % Export the data to Excel
    filename = 'data_points.xlsx';  % Excel file name
    try
        writecell(data, filename, 'Sheet', 1, 'Range', 'A1');  % Write data to Excel
        disp('Data exported to Excel successfully.');
    catch
        errordlg('Error saving data to Excel.', 'Export Error');
    end
end

% Function to manually export table to Excel via button
function exportToExcel(dataTable)
    % 获取表格数据
    data = dataTable.Data;

    % 检查是否有数据可以导出
    if isempty(data)
        errordlg('没有数据可以导出！', '导出错误');
        return;
    end

    % 获取列名
    columnNames = dataTable.ColumnName;

    % 将列名转换为列向量，以便与数据拼接
    columnNames = columnNames(:);  % 强制转换列名为列向量

    % 确保列数与数据列一致
    [numRows, numCols] = size(data);
    numColsColumnNames = numel(columnNames);  % 列名的数量

    % 如果列数不一致，填充或截断数据
    if numCols < numColsColumnNames
        % 用空单元格填充数据列
        for i = 1:numRows
            data{i, numCols+1:numColsColumnNames} = {''};
        end
    elseif numCols > numColsColumnNames
        % 如果数据列数大于列名，截断多余的列
        data = data(:, 1:numColsColumnNames);
    end

    % 合并列名与数据
    outputData = [columnNames'; data];  % 列名和数据合并为一个cell数组

    % 设置导出文件路径
    filename = 'data_points.xlsx';  % Excel 文件名
    try
        % 将列名和数据一起写入 Excel 文件
        writecell(outputData, filename, 'Sheet', 1, 'Range', 'A1');  % 将数据写入 Excel 文件
        msgbox('数据已成功导出为 Excel 文件！', '导出成功');  % 弹出成功提示框
    catch ME
        % 如果发生错误，显示详细的错误信息
        errordlg(['导出 Excel 文件时出错：' ME.message], '导出错误');  % 弹出详细错误信息
    end
end

% Function to allow user to input K() value with increment
function inputKValue(dataTable)
    % Prompt for the two numbers (K() base and initial increment value)
    prompt = {'请输入第一个数字 (固定值):', '请输入第二个数字 (起始增量值):'};
    dlgTitle = '输入 K() 值';
    numLines = 1;
    defaultAns = {'147', '900'};  % Default values for KBase and KIncrement
    answer = inputdlg(prompt, dlgTitle, numLines, defaultAns);

    if ~isempty(answer)
        % Read the input values
        KBase = str2double(answer{1});  % Fixed K() value
        KIncrement = str2double(answer{2});  % Starting increment value (e.g., 900)

        % Validate input
        if mod(KIncrement, 100) ~= 0
            errordlg('第二个数字必须是 100 的整数倍！', '输入错误');
            return;
        end

        if isnan(KBase) || isnan(KIncrement)
            errordlg('请输入有效的数字！', '输入错误');
            return;
        end

        % Save the input values for later use
        setappdata(gcf, 'KBase', KBase);  % Store KBase value
        setappdata(gcf, 'KIncrement', KIncrement);  % Store KIncrement value

        % Get the number of rows in the table
        [numRows, ~] = size(dataTable.Data);

        % For each row, update the K() value
        for i = 1:numRows
            % Calculate the K() value based on the increment logic
            increment = KIncrement + 100 * (i - 1);  % Increment by 100 for each row
            K_value = sprintf('K%d+%03d', KBase + floor(increment / 1000), mod(increment, 1000));  % Fix the formula to correctly handle carry over
            dataTable.Data{i, 1} = K_value;  % Set the K() value in the first column
        end

        msgbox('K() 值已成功输入并更新！', '输入成功');
    end
end
