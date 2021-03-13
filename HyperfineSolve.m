function [E,V,U,StateLabel,StateLabelF]=HyperfineSolve(B,KRb,L,J)


mu_B=9.27400968e-24/6.626070040e-34*1e-6*1e-4;
S=0.5;

if KRb==1,
    I=4;
    if L==0,
        A1=-285.7308;
        A2=0;
    elseif L==1 && J==0.5,
        A1=-34.523;
        A2=0;
    elseif L==1 && J==1.5,
        A1=-7.585;
        A2=-3.445;
    end;   
    gI=0.000176490;
    gJ=(J*(J+1)-S*(S+1)+L*(L+1))./(2*J*(J+1))+2*(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1));
else
    I=3/2;
    if L==0,
        Ehfs=6834.682610904;
        A1=Ehfs./(I+0.5);
        A2=0;
    elseif L==1 && J==0.5,
        A1=408.328;
        A2=0;
    elseif L==1 && J==1.5,
        A1=84.7185;
        A2=12.4965;
    end;
    gI=-0.0009951414;
    gJ=(J*(J+1)-S*(S+1)+L*(L+1))./(2*J*(J+1))+2*(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1));
end;

NI=2*I+1;
NJ=2*J+1;
Ndim=NI*NJ;
mJ=-J:J;
mI=-I:I;
StateLabel=[reshape(repmat(mI,NJ,1),Ndim,1) repmat(mJ(:),NI,1)];

F=abs(I-J):abs(I+J);
StateLabelF=[];
for nn=1:numel(F),
    mF=(-F(nn):F(nn))';
    StateLabelF=[StateLabelF;[repmat(F(nn),numel(mF),1) mF]];
end;

U=zeros(Ndim);
for a=1:Ndim,
    F=StateLabelF(a,1);
    mF=StateLabelF(a,2);
    for b=1:Ndim,
        mI=StateLabel(b,1);
        mJ=StateLabel(b,2);
        if abs(mI+mJ)>F || (mI+mJ)~=mF,
            continue;
        else
            U(a,b)=ClebschGordan(I,J,F,mI,mJ,mF);
        end;
    end;
end;


%% Bare Hamiltonian
% H=zeros(Ndim);
% F=abs(I-J):abs(I+J);
% for a=1:Ndim,
%     mI1=StateLabel(a,1);
%     mJ1=StateLabel(a,2);
%     for b=1:Ndim,
%         mI2=StateLabel(b,1);
%         mJ2=StateLabel(b,2);
%         for c=1:numel(F),
%             if abs(mI1+mJ1)>F(c) || abs(mI2+mJ2)>F(c) || (mI1+mJ1)~=(mI2+mJ2),
%                 continue;
%             end;
%                 CG=ClebschGordan(I,J,F(c),mI1,mJ1,mI1+mJ1).*ClebschGordan(I,J,F(c),mI2,mJ2,mI2+mJ2);
%                 K=F(c)*(F(c)+1)-I*(I+1)-J*(J+1);
%                 if J>=1.0 && I>=1 && L>0,
%                     tmp=CG*(A1/2*K+A2*(1.5*K*(K+1)-2*I*J*(I+1)*(J+1))./(4*I*J*(2*I-1)*(2*J-1)));
%                 else
%                     tmp=CG*A1/2*K;
%                 end;
%                 H(a,b)=H(a,b)+tmp;
%         end;
%     end;
% end;

IdotJ=0.5*(StateLabelF(:,1).*(StateLabelF(:,1)+1)-I*(I+1)-J*(J+1));
if J>=1.0 && I>=1 && L>0,
    H=A1*IdotJ+A2.*(3*IdotJ.^2+3/2.*IdotJ-I*(I+1)*J*(J+1))./(2*I*(2*I-1)*(2*J-1));
else
    H=A1*IdotJ;
end;
H=diag(H);
H=U'*H*U;



%% Zeeman field
E=zeros(size(H,1),numel(B));
V=zeros(size(H,1),size(H,2),numel(B));
for nn=1:numel(B),
    HB=diag(mu_B*B(nn)*(gJ*StateLabel(:,2)+gI*StateLabel(:,1)));
    H2=H+HB;
    [V(:,:,nn),D]=eig(H2);
    E(:,nn)=diag(D);
    [E(:,nn),idx]=sort(E(:,nn));
    V(:,:,nn)=V(:,idx,nn);
end;


end

