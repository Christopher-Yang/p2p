% editErrorBar() alters the properties of objects drawn using 
% shadedErrorBar(). 
% 
%   x: handle of plot drawn using shadedErrorBar
%   col: color of error bar
%   width: width of line

function output = editErrorBar(x, col, width)
    x.mainLine.Color = col;
    x.mainLine.MarkerFaceColor = col;
    x.mainLine.LineWidth = width;
    x.mainLine.MarkerSize = 4;
    x.patch.FaceColor = col;
    x.patch.EdgeColor = 'none';
    x.patch.FaceAlpha = 0.4;
    x.patch.HandleVisibility = 'off';
    x.edge(1).Color = 'none';
    x.edge(2).Color = 'none';
    output = x;
end