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
        G = 6.67e-11;
        tide = 1.57e-7*9.81;
        
        %%Rubidium associated "constant"
        f_Rb_groundHFS = 6834.682610904e6; %hyperfine splitting between the ground states  for Rb87
        k_Rb_groundHFS = 2*pi*6834.682610904e6/299792458;
        %wavevector tranistion from gound state to excited state  for Rb87
        k_Rb_F1_F0 = 2*pi*3.842344540713058e14/299792458;
        k_Rb_F1_F1 = 2*pi*3.842345262933378e14/299792458;
        k_Rb_F1_F2 = 2*pi*3.842346832338618e14/299792458;
        
        k_Rb_F2_F1 = 2*pi*3.842276976107269e14/299792458;
        k_Rb_F2_F2 = 2*pi*3.842278485512509e14/299792458;
        k_Rb_F2_F3 = 2*pi*3.842281152034289e14/299792458;
        %frequency transition from ground state to excited state  for Rb87
        f_Rb_F1_F0 = 3.842344540713058e14;
        f_Rb_F1_F1 = 3.842345262933378e14;
        f_Rb_F1_F2 = 3.842346832338618e14;
        
        f_Rb_F2_F1 = 3.842276976107269e14;
        f_Rb_F2_F2 = 3.842278485512509e14;
        f_Rb_F2_F3 = 3.842281152034289e14;
        
       %wavevector tranistion from gound state to excited state for Rb85
      k_Rb85_F3_F2 = 2*pi*384.2290576492837e12/299792458 ;
      k_Rb85_F3_F3 = 2*pi*384.2291210491137e12/299792458 ;
      k_Rb85_F3_F4 = 2*pi*384.2292416894837e12/299792458;
      
      k_Rb85_F2_F3 = 2*pi*384.2321567815528e12/299792458 ;
      k_Rb85_F2_F2 = 2*pi*384.2320933817228e12/299792458 ;
      k_Rb85_F2_F1 = 2*pi*384.2320640082228e12/299792458 ;
       %frequency transition from ground state to excited state  for Rb85
        f_Rb85_F3_F4 =384.2292416894837e12;
        f_Rb85_F3_F3 =384.2291210491137e12;
        f_Rb85_F3_F2 =384.2290576492837e12;
        
        f_Rb85_F2_F1 = 384.2320640082228e12;
        f_Rb85_F2_F2 =384.2320933817228e12;
        f_Rb85_F2_F3 = 384.2321567815528e12;
        %hyperfine splitting between the ground states  for Rb87
        f_Rb85_groundHFS=3035.732439060e6; 
        k_Rb85_groundHFS=2*pi*3035.732439060e6/299792458;

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
        
        function r = rejectionSample(f,inputRange,Nsamples,Pmax)
            function_inputs = strsplit(regexp(func2str(f), '(?<=^@\()[^\)]*', 'match', 'once'), ',');
            if numel(function_inputs) == 1
                if nargin < 4
                    rr = linspace(inputRange(1),inputRange(2),1e3);
                    Pmax = max(f(rr));
                end
                r = [];
                while numel(r) < Nsamples
                    rTest = (inputRange(2)-inputRange(1))*rand(Nsamples,1) + inputRange(1);
                    rr = Pmax.*rand(Nsamples,1);
                    r = [r;rTest(rr<f(rTest))];  %#ok
                end
                r = r(1:Nsamples);
            else
                num_inputs = numel(function_inputs);
                if nargin < 4
                    rr = zeros(num_inputs,1e2);
                    for nn = 1:num_inputs
                        rr(nn,:) = linspace(inputRange(nn,1),inputRange(nn,2),1e2);
                    end
                    if num_inputs == 2
                        [X,Y] = meshgrid(rr(1,:),rr(2,:));
                        f_val = f(X,Y);
                        Pmax = max(f_val(:));
                    elseif num_inputs == 3
                        [X,Y,Z] = meshgrid(rr(1,:),rr(2,:),rr(3,:));
                        f_val = f(X,Y,Z);
                        Pmax = max(f_val(:));
                    end
                end
                r = zeros(0,num_inputs);
                while size(r,1) < Nsamples
                    rTest = (inputRange(:,2) - inputRange(:,1))'.*rand(Nsamples,num_inputs) + inputRange(:,1)';
                    rr = Pmax.*rand(Nsamples,num_inputs);
                    if num_inputs == 2
                        rnew = rTest(all(rr < f(rTest(:,1),rTest(:,2)),2),:);
                        r = [r;rnew];  %#ok
                    elseif num_inputs == 3
                        r = [r;rTest(all(rr < f(rTest(:,1),rTest(:,2),rTest(:,3)),2),:)];  %#ok
                    end
                end
                r = r(1:Nsamples,:);
            end
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
                loglog(f,abs(Y)/sqrt(diff(f(1:2))),'.-');
                xlabel('Frequency [Hz]');ylabel('Noise amplitude spectral density [[x units] Hz^{-1/2}]');
            elseif strcmpi(plotType,'psd')
                loglog(f,(abs(Y)/sqrt(diff(f(1:2)))).^2,'.-');
                xlabel('Frequency [Hz]');ylabel('Noise power spectral density [[x units]^2 Hz^{-1}]');
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
            if numel(dt) > 1
                dt = abs(dt(2)-dt(1));
            end
            N = size(data,1);
            f = 1/(dt)*(0:N-1)/(N);
            f = f(1:floor(N/2));
            f = f(:);
            if nargin==3
                data = const.smooth(data,sm);
            end
            Y = fft(data,[],1);
            Y = Y(1:floor(N/2),:);
            Y(1,:) = Y(1,:)/N;
            Y(2:end,:) = Y(2:end,:)/N;
        end
        
        function y = smooth(data,sm)
            %SMOOTH Subtracts smoothed version of data from data
            %
            %   y = const.smooth(data,sm) computes a smoothed version of
            %   data using smoothing length sm and removes that from the
            %   data
            if any(size(data) == 1)
                data = data(:);
            end
            if isempty(sm)
                y = data;
                return;
            end
            y = zeros(size(data));
            for nn = 1:size(data,2)
                y(:,nn) = data(:,nn) - smooth(data(:,nn),sm);
            end
        end
        
        function y = butterworth(data,width,order)
            if nargin < 3
                order = 4;
            end
            Y = fftshift(fft(data));
            N = size(Y,1);
            f = 0.5*linspace(-1,1,N)';
            F = (1 + (f*width).^(2*order)).^-1;
            y = real(ifft(ifftshift(Y.*F)));
        end

        function y = butterworth2D(data,width,order)
            if nargin < 3
                order = 4;
            end
            Y = fftshift(fft2(data));
            kx = linspace(-0.5,0.5,size(data,2));
            ky = linspace(-0.5,0.5,size(data,1));
            [KX,KY] = meshgrid(kx,ky);
            F = 1./(1 + (width^2*(KX.^2 + KY.^2)).^order);
            y = real(ifft2(ifftshift(F.*Y)));
        end
        
        function [P,f] = calcPSD(data,dt,varargin)
            %CALCPSD Calculates the one-sided power spectral density and
            %returns it and the frequency vector
            %
            %   [P,f] = CALCPSD(DATA,DT,VARARGIN) calculates the one-sided
            %   power spectral density P and returns it and the frequency
            %   vector F based on DATA and time difference DT.  DT can also
            %   be a vector of times which is used to calculate the time
            %   difference.  VARARGIN is any valid variable argument list
            %   for CALCFFT
            [Y,f] = const.calcFFT(data,dt,varargin{:});
            Y = Y/sqrt(2);
            P = abs(Y).^2./(f(2)-f(1));
        end
        
        function [R,tau] = calcAutoCorrelation(data,dt,method,varargin)
            %CALCAUTOCORRELATION Calculates the autocorrelation of the
            %signal based on the power spectral density
            %
            %   [R,TAU] = CALCAUTOCORRELATION(DATA,DT,VARARGIN) calculates
            %   the one-sided autocorrelation function R and the vector of
            %   time delays TAU based on DATA and time difference DT.
            %   VARARGIN can be any valid variable argument list for
            %   CALCFFT
            if numel(dt) > 1
                dt = abs(dt(2)-dt(1));
            end
            if strcmpi(method,'psd')
                [P,f] = const.calcPSD(data,dt,varargin{:});
                P = [P;flipud(P(2:end,:))];
                R = numel(P)/2*(f(2) - f(1))*real(ifft(P,[],1));
                N = size(R,1);
                R = R(1:floor(N/2),:);
                tau = dt*(0:(size(R,1)-1));
            elseif strcmpi(method,'direct')
                N = size(data,1);
                if numel(varargin) == 1
                    tauStep = varargin{1};
                else
                    tauStep = max(1,floor(N/100));
                end
                tau = (0:tauStep:(floor(N/2)-1))';
                R = zeros(numel(tau),size(data,2));
                for nn = 1:numel(tau)
                    tmp1 = data(1:(end-tau(nn)),:);
                    tmp2 = data((tau(nn)+1):end,:);
                    R(nn,:) = mean(tmp1.*tmp2,1);
                end
                tau = tau.*dt;
%             elseif strcmpi(method,'overlap')
%                 N = size(data,1);
%                 if numel(varargin) == 1
%                     numSegments = varargin{1};
%                 else
%                     numSegments = max(floor(N/1e4),1);
%                 end
%                 Ns = floor(N/numSegments);
%                 data = data(1:(numSegments*Ns),:);
%                 R = zeros(Ns,size(data,2));
%                 tau = (0:(Ns-1))';
%                 for col = 1:size(data,2)
%                     idx = 
%                     tmp1 = reshape(data(
%                     
%                     for nn = 1:numel(tau)
%                         tmp1 = data(1:(end-tau(nn)),col);
%                         tmp2 = data((tau(nn)+1):end,col);
%                         R(nn,:) = mean(tmp1.*tmp2,1);
%                     end
%                     tau = tau.*dt;
%                 end
            end
        end
        
        function R = compute_tf_radii(real_trap_freqs,num_atoms,scattering_length,atom_mass)
            f = 2*pi*real_trap_freqs;
            fmean = prod(f)^(1/3);
            U0 = 4*pi*const.hbar^2*scattering_length/atom_mass;
            chemical_potential = (15*num_atoms*U0/(8*pi))^(2/5)*(atom_mass*fmean^2/2)^(3/5);

            R = sqrt(2*chemical_potential./(atom_mass*f.^2));
        end
        
        function y = sinc(x)
            if numel(x) == 1
                if x == 0
                    y = 1;
                else
                    y = sin(x)/x;
                end
            else
                y = sin(x)./x;
                y(x == 0) = 1;
            end
        end
    end
    
    
end