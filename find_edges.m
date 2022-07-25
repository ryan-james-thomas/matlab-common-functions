function [edges,edge_times] = find_edges(t,data,threshold,hysteresis)
%FIND_EDGES Finds edges in the data based on a digital implementation of a
%Schmidt trigger
%
%   [EDGES,EDGE_TIMES] = FIND_EDGES(T,DATA,THRESHOLD,HYSTERESIS) Given a
%   vector of times T with data D, a threshold value about which the data
%   needs to cross to count as an edge, and a HYSTERESIS value which the
%   data needs to cross to re-arm, returns a vector of values EDGES
%   corresponding to +1 on positive edges and -1 on negative edges which is
%   the size of DATA.  Also returns EDGE_TIMES, which is a vector
%   containing the times of edges.
armed = true;
last_edge = 0;
data = data - threshold;
edges = zeros(size(data));
edge_times = [];
mm = 1;
for nn = 2:numel(data)
    if armed && data(nn) > 0 && data(nn - 1) < 0
        edges(nn) = 1;
        last_edge = 1;
        edge_times(mm) = (data(nn)*t(nn-1) - data(nn-1)*t(nn))/(data(nn) - data(nn-1));
        mm = mm + 1;
        armed = false;
    elseif armed && data(nn) < 0 && data(nn - 1) > 0
        edges(nn) = -1;
        last_edge = -1;
        edge_times(mm) = (data(nn)*t(nn-1) - data(nn-1)*t(nn))/(data(nn) - data(nn-1));
        mm = mm + 1;
        armed = false;
    elseif (data(nn) > hysteresis && last_edge == 1) || (data(nn) < hysteresis && last_edge == -1)
        armed = true;
    end
    
end

end