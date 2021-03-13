function ar = confplot(x,m,s,varargin)
%CONFPLOT Plots a shaded region corresponding to the confidence bounds of a
%given data set.  Returns the array of area objects.

x = x(:);
m = m(:);
s = s(:);

ar = area(x,[m+s,-2*s]);
set(ar(1),'facecolor','none');
set(ar,'edgecolor','none');

if numel(varargin) == 2 && strcmpi(varargin{1},'color')
    set(ar(2),'facecolor',varargin{2});
end

