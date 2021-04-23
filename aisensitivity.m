function H = aisensitivity(tau,T,f)

w = 2*pi*f;
rabi = pi/(2*tau);
TT = T+2*tau;
G = 4*1i*rabi./(w.^2-rabi.^2).*sin(w.*TT/2).*(cos(w.*TT/2)...
    + rabi*T/2*sinc(w*T/(2*pi)));
H = w.*G;


end