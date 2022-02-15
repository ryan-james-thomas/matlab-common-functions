classdef nonlinfit < FitClass
%NONLINFIT SubClass of FitClass.  Aids in fitting non-linear functions
    properties
        %% Input data
        lower
        upper
        guess
    end
    
    properties(Access = protected)
        fo      %Fit options
        ft      %Fit type
    end
    
    methods
        %NONLINFIT Defines the FitClass object with variable length
        %arguments
        %
        %obj = NONLINFIT(X,Y,DY,EX) creates a NONLINFIT object with x and y
        %values and optional y errors DY and exclusions EX.  EX can be
        %either an array of logical values or indices to exclude
        function obj = nonlinfit(varargin)
            obj=obj@FitClass(varargin{:});
        end
        
        %BOUNDS Sets the bounds and guess for the nonlinear fit
        %
        %OBJ = BOUNDS(LOWER,UPPER,GUESS) Sets the lower and upper bounds,
        %and also sets the initial starting point or guess
        %
        %OBJ = BOUNDS() Prints the lower and upper bounds as well as the
        %guess
        function varargout = bounds(obj,lower,upper,guess)
            if nargin > 1
                obj.lower = lower;
                obj.upper = upper;
                obj.guess = guess;
                varargout{1} = obj;
            else
                if ~isempty(obj.func) && isa(obj.func,'function_handle')
                    anoninputs = strsplit(regexp(func2str(obj.func), '(?<=^@\()[^\)]*', 'match', 'once'), ',');
                    fprintf(1,'Arguments:\t');
                    for nn = 1:(numel(anoninputs)-1)
                        fprintf(1,'\t%10s',anoninputs{nn});
                    end
                    fprintf(1,'\n');
                end
                fprintf(1,'Lower bounds:');
                fprintf(1,'\t%10.4g',obj.lower);
                fprintf(1,'\n');
                fprintf(1,'Upper bounds:');
                fprintf(1,'\t%10.4g',obj.upper);
                fprintf(1,'\n');
                fprintf(1,'Guess:\t\t');
                fprintf(1,'\t%10.4g',obj.guess);
                fprintf(1,'\n');
            end
        end
        
        function self = bounds2(self,varargin)
            if mod(numel(varargin),2) ~= 0
                error('Arguments must be in name/value pairs');
            elseif isempty(self.func)
                error('A function must be supplied to use this method!');
            end
            
            anoninputs = strsplit(regexp(func2str(self.func), '(?<=^@\()[^\)]*', 'match', 'once'), ',');
            N = numel(anoninputs) - 1;
            anoninputs = anoninputs(1:N);
            
%             [self.lower,self.upper,self.guess] = deal(NaN(1,N));
            for nn = 1:numel(anoninputs)
                for mm = 1:2:numel(varargin)
                    if strcmp(varargin{mm},anoninputs{nn})
                        self.lower(nn) = varargin{mm + 1}(1);
                        self.upper(nn) = varargin{mm + 1}(2);
                        self.guess(nn) = varargin{mm + 1}(3);
                    end
                end
            end
            
        end
        
        function r = get(self,p)
            anoninputs = strsplit(regexp(func2str(self.func), '(?<=^@\()[^\)]*', 'match', 'once'), ',');
            N = numel(anoninputs) - 1;
            anoninputs = anoninputs(1:N);
            for nn = 1:N
                if strcmp(p,anoninputs{nn})
                    r = self.c(nn,:);
                end
            end
        end
        
        function print(obj)
            %PRINT Prints the fit results with the bounds and guess
            obj.bounds;
            fprintf(1,'Results:\t');
            fprintf(1,'\t%10.4g',obj.c(:,1));
            fprintf(1,'\n');
        end
        
        %MAKEFITOBJECTS Creates the internal fit objects for fitting data
        function obj = makeFitObjects(obj)
            obj.fo = fitoptions('method','NonlinearLeastSquares');
            set(obj.fo,'MaxFunEvals',1500);
            set(obj.fo,'Lower',obj.lower,'Upper',obj.upper,...
                'StartPoint',obj.guess,'Exclude',obj.ex);
            if obj.useErr
                set(obj.fo,'Weight',obj.dy.^(-2));
            end
            
            obj.ft = fittype(obj.func);
        end
        
        %FIT Performs the non-linear fit using MATLAB's fit() function
        function p = fit(obj)
            obj.makeFitObjects;
            
            [fr_,gof,output] = fit(obj.x,obj.y,obj.ft,obj.fo);
            
            cvalues = coeffvalues(fr_); %Returns coefficient values in same order as in 'fittype'
            try
                cfd_int = confint(fr_,0.68); %Returns confidence interval with column i the ith coefficient and the 2 rows the lower and upper confidence bound.  68 percent cfd bound
                err_cvalues = (cfd_int(2,:)-cfd_int(1,:))/2; %Calculate the actual errors as half the difference between the two cfd bounds
                obj.c = [cvalues(:),err_cvalues(:)]; %Combine coefficients and errors into a single variable
            catch
                obj.c = cvalues(:);
            end
            
            obj.gof.dof = gof.dfe;
            obj.gof.chi2 = gof.sse/gof.dfe;
            obj.gof.prob = 1-gammainc(gof.sse/2,gof.dfe/2);
            
            obj.res = output.residuals(:);

            p = obj.c;
        end
        
        %MONTECARLO Performs Monte Carlo analysis of the fit result
        %
        %OBJ = MONTECARLO(ITER,BOOTSTRAP) Computes the covariance matrix
        %for the resulting fit using Monte Carlo techniques using ITER
        %iterations.  If BOOTSTRAP == 0, then simulates data set with noise
        %equivalent to DY.  If BOOTSTRAP == 1, then uses the bootstrap
        %method of selecting only a subset of data before fitting
        function obj = montecarlo(obj,iter,bootstrap)
            obj.makeFitObjects;
            set(obj.fo,'StartPoint',obj.c(:,1));
            cvalues = zeros(size(obj.c,1),iter);
            for nn=1:iter
                if ~bootstrap
                    yy = obj.f(obj.x)+obj.dy.*randn(size(obj.x));
                    xx = obj.x;
                else
                    idx = randi(numel(obj.y),numel(obj.y),1);
                    yy = obj.y(idx);
                    xx = obj.x(idx);
                end
                fr_ = fit(xx,yy,obj.ft,obj.fo);
                tmp = coeffvalues(fr_);
                cvalues(:,nn) = tmp(:);
            end
            
%             obj.c(:,2) = std(cvalues,0,2);
            obj.Vcov = cov(cvalues');
            obj.Vcorr = corr(cvalues');
            
        end
        
        %DF Computes the function error
        %
        %ERR = DF(X) Computes the error in the function value for every X
        %given the previously computed covariance matrix and coefficient
        %values.
        function err = df(obj,x)
            x = x(:);          
            Ncoeffs = size(obj.c,1);
            dfdc = zeros(numel(x),Ncoeffs);
            for nn=1:Ncoeffs
                c1 = obj.c(:,1);
                c2 = obj.c(:,1);
                c1(nn) = (1-1e-6)*obj.c(nn,1);
                c2(nn) = (1+1e-6)*obj.c(nn,1);
                c1cell = num2cell(c1);
                c2cell =num2cell(c2);
                f1 = obj.func(c1cell{:},x);
                f2 = obj.func(c2cell{:},x);
                dfdc(:,nn) = (f2-f1)./(c2(nn)-c1(nn));
            end
            err = sqrt(diag(dfdc*obj.Vcov*dfdc'));
        end
        
        function s = struct(self)
            s = struct@FitClass(self);
            s.lower = self.lower;
            s.upper = self.upper;
            s.guess = self.guess;
        end
        
    
    end
    
    methods(Static)
        function s = loadobj(a)
            %LOADOBJ Converts structure into class
            s = nonlinfit(a.x,a.y,a.dy,a.ex);
            s.func = a.func;
            s.useErr = a.useErr;
            s.c = a.c;
            s.Vcov = a.Vcov;
            s.Vcorr = a.Vcorr;
            s.res = a.res;
            s.gof = a.gof;
            s.lower = a.lower;
            s.upper = a.upper;
            s.guess = a.guess;
        end
    end
    
end