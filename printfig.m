function printfig(fig,filename,varargin)

filetype = 'pdf';
fileres = 300;
printmethod = 'painters';
savematlabfig = false;
fig_size = [];

for nn=1:2:numel(varargin)
    switch varargin{nn}
        case 'type'
            filetype = varargin{nn+1};
        case 'resolution'
            fileres = varargin{nn+1};
        case 'method'
            printmethod = varargin{nn+1};
        case 'savematlabfig'
            savematlabfig = varargin{nn+1};
        case 'size'
            fig_size = varargin{nn + 1};
        otherwise
            error('Option %s unsupported',varargin{nn});
    end
end

% figure(fig);

set(fig,'units','centimeters');
pos = get(fig,'position');
if isempty(fig_size)
    fig_size = pos(3:4);
end
set(fig,'paperunits','centimeters','papersize',fig_size,'paperposition',[0,0,fig_size]);

print(fig,[filename,'.',filetype],['-d',filetype],['-r',sprintf('%d',fileres)],['-',printmethod]);
if savematlabfig
    savefig(fig,[filename,'.fig']);
end