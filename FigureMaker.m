classdef FigureMaker < handle
    properties
        fig
        units
        
        axs
        dims
        
        x0
        y0
        xw
        yw
        
        size
        font_size
    end
    
    methods
        function self = FigureMaker(varargin)
            self.units = 'centimeters';
            self.font_size = 8;
            self.dims = [1,1];
            self.set(varargin{:});
            clf;
            self.axs = axes;
            delete(self.axs);
        end
        
        function self = set(self,varargin)
            if mod(numel(varargin),2) ~= 0
                error('Arguments must be in name/value pairs!');
            else
                for nn = 1:2:numel(varargin)
                    v = varargin{nn + 1};
                    switch lower(varargin{nn})
                        case 'fig'
                            if isnumeric(v)
                                self.fig = figure(v);
                            elseif ishandle(v)
                                self.fig = v;
                            end
                        case 'dims'
                            self.dims = v;
                        case 'x0'
                            self.x0 = v;
                        case 'y0'
                            self.y0 = v;
                        case 'xw'
                            self.xw = v;
                        case 'yw'
                            self.yw = v;
                        case 'size'
                            self.set_size(v);
                        case {'font_size','fontsize'}
                            self.font_size = v;
                    end
                end
            end
        end
        
        function self = refresh_axes(self)
            
        end
        
        function self = set_size(self,size_in)
            if nargin > 1
                self.size = size_in;
            end
            
            set(self.fig,'units',self.units,'paperunits',self.units);
            self.fig.Position(3:4) = self.size;
            set(self.fig,'papersize',self.size,'paperposition',[0,0,self.size]);  
        end
        
        function self = make_all_axes(self)
            for row = 1:self.dims(1)
                for col = 1:self.dims(2)
                    self.make_axes(row,col);
                end
            end
            
        end
        
        function self = make_axes(self,row,col)
            ax = axes('position',[col*self.x0 + (col - 1)*self.xw,row*self.y0 + (row - 1)*self.yw,self.xw,self.yw]);
%             ax = axes('outerposition',[(col - 1)*self.xw,(row - 1)*self.yw,self.xw,self.yw]);
            self.axs(row,col) = ax;
        end
        
        function self = label(self,x,y,txt,varargin)
            text(x,y,txt,'units','normalized','fontweight','bold','fontsize',self.font_size,varargin{:});
        end

    end
end