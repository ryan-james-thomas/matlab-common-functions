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
                for nn=1:2:numel(varargin)
                    switch lower(varargin{nn})
                        case 'includezero'
                            includeZero = varargin{nn+1};
                        otherwise
                            error('Option %s not recognized',varargin{nn});
                    end
                end
            end
            Nplot = 1e3;
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
            xplot = xplot(:);
            yplot = obj.f(xplot);
            subplot(3,1,1:2);
            if obj.useErr
                errorbar(obj.x(~obj.ex),obj.y(~obj.ex),obj.dy(~obj.ex),'o','markersize',4);
            else
                plot(obj.x(~obj.ex),obj.y(~obj.ex),'o','markersize',4);
            end
            hold on;
            plot(xplot,yplot,'r-');
%             hold off;
            subplot(3,1,3);
            errorbar(obj.x(~obj.ex),obj.res,ones(size(obj.x(~obj.ex))),'o','markersize',4);
            hold on;
            grid on;
        end
        
        %F Returns the fit function value
        %
        %F(X) Returns the fit function value for a given X with the current
        %coefficients
        function v = f(obj,x)
            x = x(:);
            cc = num2cell(obj.c(:,1));
            v = obj.func(cc{:},x);
        end
    end
    
    methods (Abstract)
        p = fit(obj);
    end
    
    
end