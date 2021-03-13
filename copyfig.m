function copyfig(oldFig,newFig,position)

if isa(oldFig,'numeric')
    oldFig = figure(oldFig);
end

if isa(newFig,'numeric')
    newFig = figure(newFig);
end

figure(newFig);

axOld = get(oldFig,'Children');
for nn=1:numel(axOld)
    oldPosition = axOld(nn).Position;
    newPosition = oldPosition;
    newPosition([1,3]) = oldPosition([1,3])*position(3);
    newPosition([2,4]) = oldPosition([2,4])*position(4);
    newPosition(1:2) = newPosition(1:2) + position(1:2);
    if strcmpi(axOld(nn).Type,'axes')
        hNew(nn) = axes('position',newPosition);
    elseif strcmpi(axOld(nn).Type,'colorbar')
        hNew(nn) = colorbar('position',newPosition);
    else 
        continue;
    end
    copyobj(axOld(nn).Children,hNew(nn));
    p = properties(axOld(nn));
    for mm=1:numel(p)
        if includeProperties(p{mm})
            prop = get(axOld(nn),p{mm});
            if isa(prop,'handle')
                set(hNew(nn),p{mm},copy(prop));
            else
                set(hNew(nn),p{mm},prop);
            end
        end
    end
end


end

function ex = includeProperties(p)

props = {'ALim','AlimMode','BoxStyle','Box','CLim','CLimMode','Color',...
         'ColorOrder','ColorOrderIndex','DataAspectRatio','DataAspectRatioMode',...
         'FontAngle','FontName','FontSize',...
         'FontUnits','FontWeight','LineStyleOrder','LineStyleOrderIndex',...
         'LineWidth','TickDir','TickDirMode','TickLabelInterpreter',...
         'TickLength','Title','TitleFontSizeMultiplier','TitleFontWeight',...
         'Units','View','XAxisLocation','XColor','XColorMode','XDir','XGrid',...
         'XLabel','XLim','XLimMode','XMinorGrid','XMinorTick','PlotBoxAspectRatio',...
         'PlotBoxAspectRatioMode','XScale',...
         'XTick','XTickLabel','XTickLabelMode','XTickLabelRotation',...
         'XTickMode','YAxisLocation','YColor','YColorMode','YDir','YGrid',...
         'YLabel','YLim','YLimMode','YMinorGrid','YMinorTick','YScale',...
         'YTick','YTickLabel','YTickLabelMode','YTickLabelRotation',...
         'YTickMode','ZAxisLocation','ZColor','ZColorMode','ZDir','ZGrid',...
         'ZLabel','ZLim','ZLimMode','ZMinorGrid','ZMinorTick','ZScale',...
         'ZTick','ZTickLabel','ZTickLabelMode','ZTickLabelRotation',...
         'ZTickMode'};
for nn=1:numel(props)
    if strcmpi(p,props{nn})
        ex = true;
        return;
    end
end
ex = false;
         
         

end

function ex = excludeProperties(p)

props = {'Position','OuterPosition','ActivePositionProperty','Parent',...
         'Children','CurrentPoint','XAxis','YAxis','ZAxis','Selected',...
         'BeingDeleted','BusyAction','HandleVisibility','HitTest',...
         'Interruptible','UIContextMenu','Type'};
for nn=1:numel(props)
    if strcmpi(p,props{nn})
        ex = true;
        return;
    end
end
ex = false;
         
         

end