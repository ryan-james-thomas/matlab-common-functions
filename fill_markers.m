function fill_markers(ax)

ch = ax.Children;

for nn = 1:numel(ch)
    try
        set(ch(nn),'MarkerFaceColor',ch(nn).Color);
    catch
    end
end