classdef AtomInterferometer < handle
    %ATOMINTERFEROMETER Defines a class that allows for calculating various
    %things of interest regarding atom interferometers, such as pahse noise
    %and sensitivity functions
    
    methods(Static)
        function H = sensitivity(tau,T,f)
            %SENSITIVITY Calculates the weighting/sensitivity function for
            %an atom interferometer
            %
            %   H = SENSITIVITY(TAU,T,F) calculates the sensitivity for a
            %   pi/2 pulse time TAU, a pulse separation T, and frequency
            %   vector F (in real frequency in Hz)
            w = 2*pi*f;
            rabi = pi/(2*tau);
            TT = T+2*tau;
            G = 4*1i*rabi./(w.^2-rabi.^2).*sin(w.*TT/2).*(cos(w.*TT/2)...
                + rabi*T/2*sinc(w*T/(2*pi)));
            H = w.*G;
        end
        
        function noise = calcPhaseNoise(f,psd,tau,T)
            %CALCPHASENOISE Calculates the contribution to interferometer
            %phase noise from a measurement of a power spectral density
            %
            %   NOISE = CALCPHASENOISE(F,PSD,TAU,T) returns the estimated
            %   standard deviation NOISE for each pulse separation T given
            %   frequencies F and power spectral density PSD.  The PSD and
            %   F are real-frequency quantities.
            noise = zeros(numel(T),1);
            for nn = 1:numel(T)
                H = AtomInterferometer.sensitivity(tau,T(nn),f);
                idx = ~isnan(H) & ~isnan(psd) & ~isinf(H) & ~isinf(psd);
                noise(nn) = sqrt(trapz(f(idx),abs(H(idx)).^2.*psd(idx)));
            end
        end
        
    end
    
end