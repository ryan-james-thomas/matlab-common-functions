function [P3D,x,z]=InvAbel(ImageIn,roi_row,roi_col,x_center,z_center,AutoCenterEnable)
% This program computes the inverse Abel transformation of a given 2D image
% using the BASEX method.  The inputs are the input image (ImageIn), the
% region of interest in the form of roi_row and roi_col, and the x and z
% centers of the scattering pattern/cloud (in pixels).  The program returns
% the 3D scattering pattern (P3D) in (r,z) coordinates as well as the
% centered r and z coordinates.  The pattern is limited to the size of the
% region of interest.  This method of inversion produces artifacts near
% (x,r)=0 which appear in the standalone code provided by the authors of
% the BASEX method.  The program allows for asymmetry in the z-direction
    
%% Crop input image and define coordinates and other variables
ImageIn=ImageIn(roi_row,roi_col);
x_center=round(x_center);
z_center=round(z_center);

x_center=x_center-roi_row(1)+1;
z_center=z_center-roi_col(1)+1;

[Nx,Nz]=size(ImageIn);
Nk=ceil(Nx/2-5);        %Number of basis functions

z=(1:Nz)-z_center;
x=(1:Nx)-x_center;
k=linspace(0,Nk-1,Nk).';

q1=0.01;        %X inversion regularization parameter   

%% Load basis functions.
%  To improve speed the x basis functions X_k(x) are save to a .mat file.
%  The full X_k(x) are much larger than should be needed for any image.

try 
    VAR=load('InvAbelXFuncsFull');
    Xfull=VAR.Xfull;
catch err
    if (strcmp(err.identifier,'MATLAB:load:couldNotReadFile'))
        Xfull=GenXFuncs;
        save('InvAbelXFuncsFull','Xfull');
    else
        rethrow(err);
    end;
end;


%% Initial matching of image to Xfuncs using supplied center
FullCenter=(size(Xfull,2)-1)/2+1;
X=Xfull(1:Nk,(FullCenter-x_center+1):(FullCenter-x_center+Nx));   %Center and truncate basis functions
Rx=RFunc(k.^2,x);

Xinv=inv(X*(X')+q1^2*eye(Nk,Nk));
C=Xinv*X*(ImageIn);

AbelProj=(X'*C);

%% Determine new x-center by comparing with centre of mass in x direction
if ~exist('AutoCenterEnable','var'),
    AutoCenterEnable=1;
end;

xx=(1:Nx)';
% COMx=floor(sum(xx.*sum(ImageIn,2))./sum(sum(ImageIn))); %Determine centre of mass
[M,COMx]=max(sum(ImageIn,2));
COMx=floor(COMx);
if  (COMx~=x_center) && AutoCenterEnable,
    x_center=COMx;
    X=Xfull(1:Nk,(FullCenter-x_center+1):(FullCenter-x_center+Nx));
    x=(1:Nx)-x_center;
    Rx=RFunc(k.^2,x);
    Xinv=inv(X*(X')+q1^2*eye(Nk,Nk));
    C=Xinv*X*(ImageIn);
    AbelProj=(X'*C);
end;

P3D=(Rx'*C);    %Determine distribution in (r,z) space

%% QC plots

% figure(1);clf;
% subplot(2,1,1);
% plot(x,sum(ImageIn,2),'b-',x,sum(AbelProj,2),'r--');
% ylim([0,Inf]);
% subplot(2,1,2);
% plot(z,sum(ImageIn,1),'b-',z,sum(AbelProj,1),'r--');
% ylim([0,Inf])
% % plot(1:Nx,sum(ImageIn(:,:),2),'b-',1:Nx,sum(AbelProj(:,:),2),'r--');
% % ylim([0,Inf]);
% % subplot(2,1,2);
% % plot(1:Nz,sum(ImageIn,1),'b-',1:Nz,sum(AbelProj,1),'r--');
% % ylim([0,Inf]);
% 
% figure(2);clf;
% subplot(1,2,1);
% imagesc(z,x,ImageIn,[0,1]);
% subplot(1,2,2);
% imagesc(z,x,AbelProj,[0,1]);
% figure(2);clf;
% subplot(1,2,1);
% imagesc(z,x,ImageIn,[0,0.1]);
% subplot(1,2,2);
% imagesc(z,x,P3D,[0,0.1]);
% 
% figure(3);clf;
% imagesc(z,x,P3D,[0,.1]);
% 
% figure(4);clf;
% subplot(2,1,1);
% plot(x,sum(P3D,2),'b-');
% % plot(AbelImage)
% ylim([0,Inf]);
% subplot(2,1,2);
% plot(z,sum(P3D,1),'b-');
% ylim([0,Inf]);


end


function X=GenXFuncs

Nx=513; %Must be odd
Nk=150; %Run time increase with approx Nk^4

x=linspace(-(Nx-1)/2,(Nx-1)/2,Nx);
k=linspace(0,Nk-1,Nk).';
X=zeros(Nk,Nx);
for nk=1:Nk,
    R=RFunc(0:(k(nk).^2),x);
    temp=flipud(R./repmat(GammaCoeff(0:(k(nk).^2)),1,Nx));
    X(nk,:)=sqrt(pi)*GammaCoeff(k(nk).^2).*sum(repmat(AlphaCoeff(0:(k(nk).^2)),1,Nx).*temp,1);
end;


end

function Alpha=AlphaCoeff(l)

l=l(:);
N=numel(l);
Alpha=zeros(size(l));

for n=1:N,
    Alpha(n)=prod(1-1./(2*(1:l(n))));
end;

end

function G=GammaCoeff(l)

l=l(:);
G=zeros(size(l));

for n=1:numel(l),
    if l(n)<=100,
        G(n)=exp(l(n)).*l(n).^(-l(n)).*factorial(l(n));
    else
        G(n)=sqrt(2*pi*l(n));
    end
end

    
end

function R=RFunc(n,u)

N=repmat(n(:),1,numel(u))+1*eps;
U=repmat(u(:).',numel(n),1)+1*eps;

% nlog=n.*log(n);
% nlog(n==0)=0;
% Nlog=repmat(nlog(:),1,numel(u));
% ulog=log(u.^2);
% ulog(u==0)=0;
% Ulog=repmat(ulog,numel(n),1);

% R=exp(N-U.^2+N.*Ulog-Nlog);
R=exp(N-U.^2+N.*log(U.^2)-N.*log(N));

end




