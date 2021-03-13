function [Freq,FreqLin]=BreitRabi(B,KRb,F,mF)
%Calculates Rb and K ground state hyperfine energies in MHz according to
%the Breit-Rabi formula.  Fundamental constants are 2012 NIST values, and
%atomic constants are from Steck (Rb) and Tiecke (K)
if nargin<4,
    mF=-F:F;
end;

mu_B=9.27400968e-28/6.62606957e-34*1e-6;                %MHz/G
gS=2.0023193043622;
mElectron=9.10938291e-31;   %kg

if KRb==1,
    mK=39.96399848*1.660538921e-27;   %kg
    gL=1-mElectron./mK;    
    I=4;
    Ehfs=-285.7308*(I+0.5);
    L=0;J=0.5;S=0.5;
    gI=0.000176490;
    gJ=gL.*(J*(J+1)-S*(S+1)+L*(L+1))./(2*J*(J+1))+gS*(J*(J+1)+S*(S+1)-L*(L+1))./(2*J*(J+1));
    gF=gJ.*(F*(F+1)-I*(I+1)+J*(J+1))./(2*F.*(F+1))+gI.*(F*(F+1)+I*(I+1)-J*(J+1))./(2*F.*(F+1));
else
    mRb=86.909180520*1.660538921e-27;    %kg
    gL=1-mElectron./mRb;
    Ehfs=6834.682610904;
    I=3/2;
    L=0;J=0.5;S=0.5;
    gI=-0.0009951414;
    gJ=gL.*(J*(J+1)-S*(S+1)+L*(L+1))./(2*J*(J+1))+gS*(J*(J+1)+S*(S+1)-L*(L+1))./(2*J*(J+1));
    gF=gJ.*(F*(F+1)-I*(I+1)+J*(J+1))./(2*F.*(F+1))+gI.*(F*(F+1)+I*(I+1)-J*(J+1))./(2*F.*(F+1));
end;

Freq=zeros(numel(B),numel(mF));
FreqLin=zeros(numel(B),numel(mF));
for n=1:numel(mF),

    x=(gJ-gI).*mu_B.*B./Ehfs;

    if KRb==1,
        if F==9/2,
            if mF(n)==-9/2,
                Freq(:,n)=Ehfs*(I/(2*I+1))-0.5*(gJ+2*I*gI).*mu_B.*B;
            elseif mF(n)==9/2,
                Freq(:,n)=Ehfs*(I/(2*I+1))+0.5*(gJ+2*I*gI).*mu_B.*B;
            else
                Freq(:,n)=-Ehfs/(2*(2*I+1))+gI*mu_B*mF(n)*B+Ehfs/2.*sqrt(1+4*mF(n).*x./(2*I+1)+x.^2);
            end;
%             Freq(:,n)=Freq(:,n)-Ehfs/2*(1-1./(2*I+1));
%             Freq(:,n)=-Freq(:,n);
            FreqLin(:,n)=mF(n)*gF*mu_B*B+Ehfs/2.*(1-1./(2*I+1));
        elseif F==7/2,
            Freq(:,n)=-Ehfs/(2*(2*I+1))+gI*mu_B*mF(n)*B-Ehfs/2.*sqrt(1+4*mF(n).*x./(2*I+1)+x.^2);
%             Freq(:,n)=Freq(:,n)-Ehfs/2*(-1-1./(2*I+1));
%             Freq(:,n)=-Freq(:,n);
            FreqLin(:,n)=mF(n)*gF*mu_B*B+Ehfs/2.*(-1-1./(2*I+1));
        end;
    else
        if F==2,
            if mF(n)==-2,
                Freq(:,n)=Ehfs*(I/(2*I+1))-0.5*(gJ+2*I*gI).*mu_B.*B;
            elseif mF(n)==2,
                Freq(:,n)=Ehfs*(I/(2*I+1))+0.5*(gJ+2*I*gI).*mu_B.*B;
            else
                Freq(:,n)=-Ehfs/(2*(2*I+1))+gI*mu_B*mF(n)*B+Ehfs/2.*sqrt(1+4*mF(n).*x./(2*I+1)+x.^2);
            end;
%             Freq(:,n)=Freq(:,n)-Ehfs/2*(1-1./(2*I+1));
%             Freq(:,n)=-Freq(:,n);
        elseif F==1,
            Freq(:,n)=-Ehfs/(2*(2*I+1))+gI*mu_B*mF(n)*B-Ehfs/2.*sqrt(1+4*mF(n).*x./(2*I+1)+x.^2);
%             Freq(:,n)=Freq(:,n)-Ehfs/2*(-1-1./(2*I+1));
%             Freq(:,n)=-Freq(:,n);
        end; 
    end;
    
end;


