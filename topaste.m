function topaste(FigNum,Style)
% Copies the current figure as a cropped screenshot to the clipboard
% Usage: topaste('enhanced')

if exist('FigNum','var') == 1
    figure(FigNum)
end
if exist('Style','var') == 0
    Style = '';
end




% Set the figure background color to white
set(gcf, 'Color', 'white');

% Capture the screenshot
img = getframe(gcf);
screenshot = img.cdata;

% Convert the screenshot to RGB if it is grayscale
if size(screenshot, 3) == 1
    screenshot = repmat(screenshot, [1 1 3]);
end

% Check if the 'enhanced' argument is provided
if strcmpi(Style, 'enhanced')
    % Set the font size and interpreter for the axis titles
    titleFontSize = 22;
    titleInterpreter = 'latex';
    
    % Set the font size and interpreter for the axis labels and tick labels
    axisFontSize = 22;
    axisInterpreter = 'latex';

    % Set the font size and interpreter for axis labels and tick labels
    ax = gca;
    set(ax, 'FontSize', axisFontSize, 'TickLabelInterpreter', axisInterpreter);

    % Set the font size and interpreter for axis titles
    set(ax.XLabel, 'FontSize', titleFontSize, 'Interpreter', titleInterpreter);
    set(ax.YLabel, 'FontSize', titleFontSize, 'Interpreter', titleInterpreter);
    set(ax.ZLabel, 'FontSize', titleFontSize, 'Interpreter', titleInterpreter);
end

% Find the bounding box of the figure
bbox = find_figure_bbox(screenshot);

% Crop the screenshot to the bounding box
screenshot = imcrop(screenshot, bbox);

% Copy the cropped screenshot to the clipboard using imclipboard
imclipboard('copy', screenshot);

% Notify the user that the screenshot is copied to the clipboard
fprintf('Copied screenshot to clipboard.\n');
end

function bbox = find_figure_bbox(screenshot)
% Find the bounding box of the figure within the screenshot

% Convert the screenshot to grayscale
screenshot_gray = rgb2gray(screenshot);

% Find the rows and columns that contain non-white pixels
[row, col] = find(screenshot_gray < 255);

% Calculate the bounding box coordinates
xmin = min(col);
xmax = max(col);
ymin = min(row);
ymax = max(row);

% Create the bounding box
bbox = [xmin ymin xmax-xmin+1 ymax-ymin+1];
end
