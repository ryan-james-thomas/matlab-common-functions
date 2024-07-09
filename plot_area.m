function [ha,h,h2] = plot_area(x,y,dy,args)

arguments
    x
    y
    dy
    args.Color = 'b';
    args.Alpha = 1;
    args.MeanColor = 'k';
    args.PlotMean = 0
    args.PlotBounds = 0;
    args.BoundColor = 'k';
end

if size(dy,2) == 2
    dy = abs(diff(dy,1,2));
end
ha = area(x,[y + dy,-2*dy]);
set(ha(1),'facecolor','none');
set(ha(2),'facecolor',args.Color,'facealpha',args.Alpha);
set(ha,'edgecolor','none');

if args.PlotMean
    hold on
    h = plot(x,y,'-','color',args.MeanColor);
    hold off;
else
    h = [];
end

if args.PlotBounds
    hold on
    h2 = plot(x,y + [-dy,dy],'-','color',args.BoundColor);
    hold off;
else
    h2 = [];
end
