classdef linfit < FitClass
%LINFIT A SubClass of FitClass for fitting linear functions

    properties
        non_negative
    end

    methods
        %LINFIT Defines the FitClass object with variable length
        %arguments
        %
        %obj = LINFIT(X,Y,DY,EX) creates a LINFIT object with x and y
        %values and optional y errors DY and exclusions EX.  EX can be
        %either an array of logical values or indices to exclude
        function obj = linfit(varargin)
            obj=obj@FitClass(varargin{:});
            obj.non_negative = false;
        end

        %SETFITFUNC Sets the fitting function
        %
        %SETFITFUNC(FUNC) Sets the fitting function to FUNC.  FUNC must be
        %both a function handle and return a matrix of values with
        %different linear terms in each column.  So a simple line would be
        %@(x) [ones(size(x(:))) x(:)]
        function setFitFunc(obj,varargin)
            if numel(varargin) == 2 && strcmpi(varargin{1},'poly')
                if numel(varargin{2}) == 1
                    idx = 0:varargin{2};
                else
                    idx = varargin{2};
                end
                strc = cell(numel(idx),1);
                for nn = 1:numel(idx)
                    if idx(nn) == 0
                        strc{nn} = 'ones(size(x(:)))';
                    else
                        strc{nn} = sprintf('x(:).^%d',round(idx(nn)));
                    end
                end
                str = '[';
                for nn = 1:(numel(strc)-1)
                    str = [str,strc{nn},', ']; %#ok<*AGROW>
                end
                str = ['@(x) ',str,strc{end},']'];
                obj.func = str2func(str);
%                 error('Currently not supported')
            elseif numel(varargin) == 1
                func = varargin{1};
                if ~isa(func,'function_handle')
                    error('func must be a function handle');
                elseif size(func(0),2) < 1
                    error('func must be a row vector');
                end
                obj.func = func;
            end
        end
        
        %FIT Performs a linear fit using matrices
        function p = fit(obj)
            % V = diag(obj.dy(~obj.ex).^(-2));
            V = obj.dy(~obj.ex).^(-2);
            C = obj.func(obj.x(~obj.ex));
            CV = (V.*C)';
            A = CV*C;
            if obj.non_negative
                p = lsqnonneg(C,obj.y(~obj.ex));
            else
                b = CV*obj.y(~obj.ex);
                p = A\b;
            end
            p(:,1) = p(:);
            obj.Vcov = inv(A);
            p(:,2) = sqrt(diag(obj.Vcov));
            obj.Vcorr = obj.Vcov./(p(:,2)*p(:,2).');
            obj.c = p;
            
            obj.res = (obj.y - obj.f(obj.x))./obj.dy;
            obj.gof.dof = numel(obj.y(~obj.ex)) - numel(obj.func(0));
            obj.gof.chi2 = sum(obj.res(~obj.ex).^2)/obj.gof.dof;
            obj.gof.prob = 1 - gammainc(obj.gof.chi2/2,obj.gof.dof/2);
        end
        
        %F Returns the function value
        %
        %V = F(X) Returns the function value V at X
        function v = f(obj,x)
            x = x(:);
            v = obj.func(x)*obj.c(:,1);
        end
        
        %DF Computes the function error
        %
        %ERR = DF(X) Computes the error in the function value for every X
        %given the previously computed covariance matrix and coefficient
        %values.
        function err = df(obj,x)
            x = x(:);
            err = sqrt(diag(obj.func(x)*obj.Vcov*obj.func(x)'));
        end

        function s = struct(self)
            s = struct@FitClass(self);
            s.non_negative = self.non_negative;
        end
        
    end

    methods(Static)
        function s = loadobj(a)
            %LOADOBJ Converts structure into class
            s = linfit(a.x,a.y,a.dy,a.ex);
            s.func = a.func;
            s.useErr = a.useErr;
            s.c = a.c;
            s.Vcov = a.Vcov;
            s.Vcorr = a.Vcorr;
            s.res = a.res;
            s.gof = a.gof;
            s.non_negative = a.non_negative;
        end
    end

end