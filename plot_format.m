function plot_format(x,y,t,fsize,varargin)
xlabel(x,'FontSize',fsize,varargin{:});
ylabel(y,'FontSize',fsize,varargin{:});
title(t,'FontSize',fsize,varargin{:});
set(gca,'FontSize',fsize);
