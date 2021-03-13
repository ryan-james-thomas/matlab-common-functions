classdef allandev
    
    methods(Static)
        function [N,tau,y] = prep(y,varargin)
            [rows,cols] = size(y);
            if rows == 1 || cols == 1
                y = y(:);
                rows = numel(y);
            end
            
            tau = [];
            if nargin >= 2
                if numel(varargin{1}) == 1
                    N = varargin{1};
                    if isempty(N)
                        N = min(100,rows);
                    else
                        N = min(varargin{1},rows);
                    end
                else
                    tau = varargin{1};
                    N = numel(tau);
                end
            else
                N = min(100,rows);
            end
            
            if isempty(tau)
                tau = ceil(logspace(0,log10(rows/2),N));
            end
            
            
        end
        
        function [adev,tau] = normal(y,varargin)
%             [N,tau,y] = allandev.prep(y,varargin{:});
            [rows,cols] = size(y);
            tau = unique(floor(rows./(rows:-1:1)));
            N = numel(tau);
            
            avar = zeros(N,cols);
            for cc = 1:cols
                for nn = 1:N
                    B = floor(rows/tau(nn));
                    s = reshape(y(1:(B*tau(nn)),cc),tau(nn),B);
                    avar(nn,cc) = 0.5*var(diff(mean(s,1)));
                end
            end
            
            adev = sqrt(avar);
            
        end
        
        function [adev,tau] = overlap(y,varargin)
            [Ntau,tau,y] = allandev.prep(y,varargin{:});
            [~,cols] = size(y);
            x = cumsum(y,1);
            x = [zeros(1,cols);x];
            N = size(x,1);
            
            avar = zeros(Ntau,cols);
            for cc = 1:cols
                for nn = 1:Ntau
                    i0 = 0:(N-2*tau(nn)-1);
                    i1 = tau(nn):(N-tau(nn)-1);
                    i2 = (2*tau(nn)):(N-1);
                    s0 = x(i0+1,cc);
                    s1 = x(i1+1,cc);
                    s2 = x(i2+1,cc);
                    avar(nn,cc) = 1./(2*(N-2*tau(nn)).*tau(nn).^2).*sum((s2-2*s1+s0).^2);
                end
            end
            adev = sqrt(avar);
        end
        
        function [acov,tau] = covariance(y1,y2,varargin)
            [Ntau,tau,y1] = allandev.prep(y1,varargin{:});
            [~,~,y2] = allandev.prep(y2,varargin{:});
            [~,cols] = size(y1);
            x1 = cumsum(y1,1);
            x1 = [zeros(1,cols);x1];
            x2 = cumsum(y2,1);
            x2 = [zeros(1,cols);x2];
            N = size(x1,1);
            
            acov = zeros(Ntau,cols);
            for cc = 1:cols
                for nn = 1:Ntau
                    i0 = 0:(N-2*tau(nn)-1);
                    i1 = tau(nn):(N-tau(nn)-1);
                    i2 = (2*tau(nn)):(N-1);
                    s0 = x1(i0+1,cc);
                    s1 = x1(i1+1,cc);
                    s2 = x1(i2+1,cc);
                    r0 = x2(i0+1,cc);
                    r1 = x2(i1+1,cc);
                    r2 = x2(i2+1,cc);
                    acov(nn,cc) = 1./(2*(N-2*tau(nn)).*tau(nn).^2).*sum((s2-2*s1+s0).*(r2-2*r1+r0));
                end
            end
        end
        
        function [acorr,tau] = correlation(y1,y2,varargin)
            adev1 = allandev.overlap(y1,varargin{:});
            adev2 = allandev.overlap(y2,varargin{:});
            [acov,tau] = allandev.covariance(y1,y2,varargin{:});
            acorr = acov./(adev1.*adev2);
        end
    end

end