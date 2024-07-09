classdef AtomInterferometer < handle
    %ATOMINTERFEROMETER Defines a class that allows for calculating various
    %things of interest regarding atom interferometers, such as pahse noise
    %and sensitivity functions
    
    methods(Static)
        function [H,G] = sensitivity(tau,T,f,ai_type)
            %SENSITIVITY Calculates the weighting/sensitivity function for
            %an atom interferometer
            %
            %   H = SENSITIVITY(TAU,T,F) calculates the sensitivity for a
            %   pi/2 pulse time TAU, a pulse separation T, and frequency
            %   vector F (in real frequency in Hz)
            %
            %   H = SENSITIVITY(__,AI_TYPE) calculates the sensitivity
            %   function for either 'MZ' or 'ramsey' interferometers
            w = 2*pi*f;
            rabi = pi/(2*tau);
            if nargin < 4
                ai_type = 'mz';
            end
            if strcmpi(ai_type,'mz')
                TT = T+2*tau;
                G = 4*1i*rabi./(w.^2-rabi.^2).*sin(w.*TT/2).*(cos(w.*TT/2)...
                    + rabi*T/2*sinc(w*T/(2*pi)));
            elseif strcmpi(ai_type,'ramsey')
                g1 = exp(-1i*w*tau/2).*(rabi + 1i*w.*exp(1i*w*tau))./(w.^2 - rabi.^2);
                g2 = (exp(1i*w*(T - tau/2)) - exp(1i*w*tau/2))./(1i*w);
                g2(isnan(g2)) = T - tau;
                g3 = -exp(-1i*w*(3*tau/2 - T)).*(rabi + 1i*w.*exp(1i*w*tau))./(w.^2 - rabi.^2);
                G = g1 + g2 + g3;
            else
                error('Only allowed interferometer types are ''mz'' and ''ramsey''!');
            end
            H = 1i*w.*G;
        end
        
        function noise = calcPhaseNoise(f,psd,tau,T,weights,ai_type)
            %CALCPHASENOISE Calculates the contribution to interferometer
            %phase noise from a measurement of a one-sided power spectral
            %density
            %
            %   NOISE = CALCPHASENOISE(F,PSD,TAU,T) returns the estimated
            %   standard deviation NOISE for each pulse separation T given
            %   frequencies F and one-sided power spectral density PSD.
            %   The PSD and F are real-frequency quantities.
            %
            %   NOISE = CALCPHASENOISE(__,WEIGHTs) Applies weights WEIGHTS
            %   to the PSD for each time T.  WEIGHTS should have the same
            %   length as T
            
            if nargin < 5
                weights = ones(size(T));
                ai_type = 'mz';
            elseif nargin < 6
                if isempty(weights)
                    weights = ones(size(T));
                end
                ai_type = 'mz';
            end
            noise = zeros(numel(T),1);
            for nn = 1:numel(T)
                H = AtomInterferometer.sensitivity(tau,T(nn),f,ai_type);
                idx = ~isnan(H) & ~isnan(psd) & ~isinf(H) & ~isinf(psd);
                noise(nn) = sqrt(trapz(f(idx),abs(H(idx)).^2.*psd(idx).*weights(nn)));
            end
        end

        function [H,G] = ramsey_sensitivity(tau,T,f)
            w = 2*pi*f;
            rabi = pi/(2*tau);
            g1 = exp(-1i*w*tau/2).*(rabi + 1i*w.*exp(1i*w*tau))./(w.^2 - rabi.^2);
            g2 = (exp(1i*w*(T - tau/2)) - exp(1i*w*tau/2))./(1i*w);
            g2(isnan(g2)) = T - tau;
            g3 = -exp(-1i*w*(3*tau/2 - T)).*(rabi + 1i*w.*exp(1i*w*tau))./(w.^2 - rabi.^2);
            G = g1 + g2 + g3;
            H = 1i*w.*G;
        end
        
    end
    
end