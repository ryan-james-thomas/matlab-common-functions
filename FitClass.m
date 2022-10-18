classdef FitClass < handle
%FITCLASS Defines the abstract class FitClass which can be used for linear
%and non-linear fitting with appropriate extensions of the class
    properties
        %% Input data
        x       %x-data
        y       %y-data
        dy      %y-err
        ex      %excluded points
        func    %the fit function 
        useErr  %Use the supplied errors
        %% Fitting data
        c       %coefficients
        Vcov    %covariance matrix
        Vcorr   %correlation matrix
        res     %normalized residuals
        gof     %Goodness of fit structure
    end
    
    
    methods
        %FITCLASS Defines the FitClass object with variable length
        %arguments
        %
        %obj = FITCLASS(X,Y,DY,EX) creates a FitClass object with x and y
        %values and optional y errors DY and exclusions EX.  EX can be
        %either an array of logical values or indices to exclude
        function obj = FitClass(varargin)
            obj.set(varargin{:});
        end
        
        %SET Sets the fit data
        %
        %obj = SET(X,Y,DY,EX) sets the fit data for a FITCLASS object with x and y
        %values and optional y errors DY and exclusions EX.  EX can be
        %either an array of logical values or indices to exclude
        function obj = set(obj,x,y,dy,ex)
            if nargin == 3
                obj.x = x(:);
                obj.y = y(:);
                obj.dy = 1e-2*range(obj.y).*ones(size(obj.x));
                obj.useErr = false;
                obj.ex = false(size(obj.x));
            elseif nargin >= 4
                obj.x = x(:);
                obj.y = y(:);
                obj.dy = dy(:);
                obj.useErr = true;
                if nargin == 4
                    obj.ex = false(size(obj.x));
                elseif isa(ex,'logical')
                    obj.ex = ex(:);
                else
                    obj.ex = false(size(obj.x));
                    obj.ex(ex) = true;
                end
            end
            obj.ex = obj.ex | isnan(obj.y) | isinf(obj.y);
            obj.checkSizes();
        end
        
        %CHECKSIZES Checks the sizes of the x, y, dy, and ex data to make
        %sure that they are all the same size
        function obj = checkSizes(obj)          
            if numel(obj.x) ~= numel(obj.y)
                error('Number of x points must equal number of y points');
            elseif numel(obj.x) ~= numel(obj.dy)
                error('Number of x points must equal number of dy points');
            elseif numel(obj.x) ~= numel(obj.ex)
                error('Number of x points must equal number of ex points');
            end
        end
        
        %SET.X Sets the x data
        function set.x(obj,x)
            obj.x = x(:);
        end
        
        %SET.Y Sets the y data
        function set.y(obj,y)
            obj.y = y(:);
        end
        
        %SET.DY Sets the dy data, or y error.  Can either be an array with
        %the same number of elements as y, or a single value which is then
        %replicated
        function set.dy(obj,dy)
            if numel(dy) == 1
                obj.dy = dy*ones(size(obj.y(:))); %#ok
            else
                obj.dy = dy(:);
            end
            obj.useErr = true;  %#ok
        end
        
        %SET.EX Sets the exclusion data.  Can either be an array of logical
        %values or a set of indices
        function set.ex(obj,ex)
            if isa(ex,'logical')
                obj.ex = ex(:);
            else
                obj.ex = false(size(obj.x));    %#ok
                obj.ex(ex) = true;
            end
            obj.checkSizes();
        end
        
        %SETFITFUNC Sets the fitting function as an anonymous function.
        %This function should be of the form @(a,b,c,...,x) where a,b,c,...
        %are the coefficients and x is the x data
        function setFitFunc(obj,func)
            if ~isa(func,'function_handle')
                error('func must be a function handle');
            end
            obj.func = func;
        end

        %PLOT Plots the data and fit function together along with residuals
        function plot(obj,varargin)
            if mod(numel(varargin),2)~=0
                error('Variable arguments must occur in name/value pairs!');
            else
                includeZero = false;
                plotResiduals = true;
                for nn=1:2:numel(varargin)
                    switch lower(varargin{nn})
                        case 'includezero'
                            includeZero = varargin{nn+1};
                        case 'plotresiduals'
                            plotResiduals = varargin{nn+1};
                        otherwise
                            error('Option %s not recognized',varargin{nn});
                    end
                end
            end
            Nplot = 1e3;
            if Nplot > numel(obj.x)
                xExt = 0.1*range(obj.x(~obj.ex));
                xMin = min(obj.x(~obj.ex));
                xMax = max(obj.x(~obj.ex));
                if includeZero
                    if xMin > 0
                        xplot = linspace(0,xMax+xExt,Nplot);
                    elseif xMax < 0
                        xplot = linspace(xMin-xExt,0,Nplot);
                    else
                        xplot = linspace(xMin-xExt,xMax+xExt,Nplot);
                    end
                else
                    xplot = linspace(xMin-xExt,xMax+xExt,Nplot);
                end
            else
                xplot = obj.x(~obj.ex);
            end
            xplot = xplot(:);
            yplot = obj.f(xplot);
            if plotResiduals
                subplot(3,1,1:2);
            end
            if obj.useErr
                errorbar(obj.x(~obj.ex),obj.y(~obj.ex),obj.dy(~obj.ex),'o','markersize',4);
            else
                plot(obj.x(~obj.ex),obj.y(~obj.ex),'o','markersize',4);
            end
            hold on;
            plot(xplot,yplot,'r-');
%             hold off;
            if plotResiduals
                subplot(3,1,3);
                errorbar(obj.x(~obj.ex),obj.res,ones(size(obj.x(~obj.ex))),'o','markersize',4);
                hold on;
                grid on;
            end
        end
        
        %F Returns the fit function value
        %
        %F(X) Returns the fit function value for a given X with the current
        %coefficients
        function v = f(obj,x,cc)
            x = x(:);
            if nargin == 2
                cc = num2cell(obj.c(:,1));
                v = obj.func(cc{:},x);
            else
                cc = num2cell(cc);
                v = obj.func(cc{:},x);
            end
        end
        
        function self = montecarlo(self,iter,bootstrap)
            %MONTECARLO Performs Monte Carlo analysis of the fit result
            %
            %   OBJ = MONTECARLO(ITER,BOOTSTRAP) Computes the covariance
            %   matrix for the resulting fit using Monte Carlo techniques
            %   using ITER iterations.  If BOOTSTRAP == 0, then simulates
            %   data set with noise equivalent to DY.  If BOOTSTRAP == 1,
            %   then uses the bootstrap method of selecting only a subset
            %   of data before fitting
            
            s = struct(self);
            coeffs = zeros(iter,size(s.c,1));
            
            for nn = 1:iter
                if nargin < 3 || bootstrap
                    idx = randi(numel(self.y),numel(self.y),1);
                    self.set(s.x(idx),s.y(idx),s.dy(idx),s.ex(idx));
                else
                    self.dy = self.f(self.x) + self.dy*randn(size(self.dy));
                end
                tmp = self.fit;
                coeffs(nn,:) = tmp(:,1)';
            end
            
            self.Vcov = cov(coeffs);
            self.Vcorr = corr(coeffs);
            self.gof = s.gof;
            self.c = s.c;
            self.res = s.res;
            self.x = s.x;
            self.y = s.y;
            self.dy = s.dy;
            self.ex = s.ex;
        end
        
        function s = struct(self)
            %STRUCT Creates a struct representing the object
            
            s.x = self.x;
            s.y = self.y;
            s.dy = self.dy;
            s.ex = self.ex;
            s.func = self.func;
            s.useErr = self.useErr;
            s.c = self.c;
            s.Vcov = self.Vcov;
            s.Vcorr = self.Vcorr;
            s.res = self.res;
            s.gof = self.gof;
        end
        
        function s = saveobj(self)
            %SAVEOBJ Creates a structure representing the object for saving
            s = struct(self);
        end
    end
    
    methods (Abstract)
        p = fit(obj);
    end
    
    methods(Static)
        function s = loadobj(a)
            %LOADOBJ Converts structure into class
            s = FitClass(a.x,a.y,a.dy,a.ex);
            s.func = a.func;
            s.useErr = a.useErr;
            s.c = a.c;
            s.Vcov = a.Vcov;
            s.Vcorr = a.Vcorr;
            s.res = a.res;
            s.gof = a.gof;
        end
    end
    
    
end
