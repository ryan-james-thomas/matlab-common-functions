classdef const < handle
    properties(Constant)
        %% Universal constants
        amu=1.660538921e-27;
        mK = 39.963998166*const.amu;  %Mass of 40K
        mRb = 86.909180527*const.amu;    %Mass of 87Rb
        mK41=40.96182576*const.amu;   %Mass of 41K
        mRb85 = 84.911789738*const.amu;
        mK39 = 38.96370668*const.amu;

        e=1.609e-19;    %C
        h=6.62606957e-34;
        hbar=const.h/(2*pi);         %Js
        kb=1.3806488e-23;      %J/K
        c=299792458;           %m/s
        muB=9.27400968e-24;   %J/T
        muBh=const.muB*1e-4*1e-6/const.h; %muB/h in MHz/G
        aBohr=0.52917721092e-10;    %m
        mu0=4*pi*1e-7;         %N/A^2
        eps0=1./(const.mu0*const.c^2);
        gS=2.0023193043622;
        alpha=1/137.035999074;
        EHartree=const.hbar*const.c*const.alpha/const.aBohr;
        
        g=9.81;
    end
    
    methods(Static)
        function a=cm2K
            a=1e2*const.h*const.c/const.kb;
        end
        
        function [x,y,dy] = randStat(xin,yin,N)
            xin = xin(:);
            yin = yin(:);
            [x,idx] = sort(xin);
            x = x(1:N:end);
            yBlock = reshape(yin(idx),N,round(numel(xin)/4))';
            y = mean(yBlock,2);
            dy = std(yBlock,0,2);
        end
        
        function [randOut,reverseKey] = randomize(in)
            in = in(:);
            k = randperm(numel(in));
            randOut = in(k);
            [~,reverseKey] = sort(k);
        end
        
        function x = rejectionSample(f,inputRange,Nsamples)
            xx = linspace(inputRange(1),inputRange(2),1e3);
            Pmax = max(f(xx));
            x = [];
            while numel(x)<Nsamples
                xTest = (inputRange(2)-inputRange(1))*rand(Nsamples,1)-inputRange(1);
                r = Pmax.*rand(Nsamples,1);
                x = [x;xTest(r<f(xTest))];  %#ok
            end
            x = x(1:Nsamples);
        end
        
        function plotfft(data,dt,sm,plotType)
            %PLOTFFT Plots FFT of the data
            %
            %   const.plotfft(data,dt) plots the FFT of the data with time
            %   step dt.  RMS power is plotted
            %   const.plotfft(data,dt,sm) plots the FFT of the data with
            %   time step dt with the smoothed signal removed.  Useful for
            %   removing DC offsets and other low frequency components.
            %   RMS power is plotted
            %   const.plotfft(data,dt,sm,plotType) plots the FFT of the
            %   data with different y axes.  Allowed values are 'pow' for
            %   RMS power, 'amp' for RMS amplitude, and 'nsd' for noise
            %   spectral density.
            if nargin>=3
                [Y,f]=const.calcFFT(data,dt,sm);
            else
                [Y,f]=const.calcFFT(data,dt);
            end
            Y=Y/sqrt(2);    %Appropriate for noise calculations as RMS/average values
            if nargin<4 || strcmpi(plotType,'pow')
                plot(f,abs(Y).^2,'.-');
                xlabel('Frequency [Hz]');ylabel('Average Power [arb. units]');
            elseif strcmpi(plotType,'amp')
                plot(f,abs(Y),'.-');
                xlabel('Frequency [Hz]');ylabel('RMS Amplitude [V]');
            elseif strcmpi(plotType,'nsd')
                plot(f,abs(Y)/sqrt(diff(f(1:2))),'.-');
                xlabel('Frequency [Hz]');ylabel('Noise spectral density [VHz^{-1/2}]');
            end
            xlim([0,max(f)]);
        end
        
        function [Y,f]=calcFFT(data,dt,sm)
            %CALCFFT Calculates the FFT of the data
            %
            %   const.calcFFT(data,dt) calculates the FFT of the data with 
            %   time step dt.
            %
            %   const.calcFFT(data,dt,sm) calculates the FFT of the data 
            %   with time step dt with the smoothed signal removed.  
            N=size(data,1);
            f=1/(dt)*(0:N-1)/(N);
            f=f(1:floor(N/2));
            f=f(:);
            if nargin==3
                data=servo.smooth(data,sm);
            end
            Y=fft(data,[],1);
            Y=Y(1:floor(N/2),:);
            Y(1,:) = Y(1,:)/N;
            Y(2:end,:) = 2*Y(2:end,:)/N;
        end
    end
    
    
end