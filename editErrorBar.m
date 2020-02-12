function output = editErrorBar(x, col, width)
% changes the plotting properties of shadedErrorBar.m without changing any
% of the licensed code
    x.mainLine.Color = col;
    x.mainLine.MarkerFaceColor = col;
    x.mainLine.LineWidth = width;
    x.mainLine.MarkerSize = 4;
    x.patch.FaceColor = col;
    x.patch.EdgeColor = 'none';
    x.patch.FaceAlpha = 0.15;
    x.patch.HandleVisibility = 'off';
    output = x;
end