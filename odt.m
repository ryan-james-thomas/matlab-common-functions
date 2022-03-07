classdef odt < handle
    properties
        atomType
        mass
        wavelength
        a
        useGravity
        wx
        wy
        P
    end
    
    properties(Constant)
        Rb87 = [struct('wavelength',780e-9,'decay',2*pi*6.065e6,'J',3/2);
                struct('wavelength',795e-9,'decay',2*pi*5.74e6,'J',1/2)];
%                 struct('wavelength',420e-9,'decay',2*pi*1.77e6);
%                 struct('wavelength',421.5e-9,'decay',2*pi*1.5e6);];
            
        K40 = [struct('wavelength',767e-9,'decay',2*pi*6.035e6,'J',3/2);
               struct('wavelength',770e-9,'decay',2*pi*5.956e6,'J',1/2)];

    end
    
    methods
        function obj = odt(wavelength,atomType,useGravity)
            obj.atomType = atomType;
            obj.wavelength = wavelength;
            if strcmpi(atomType,'Rb87')
                obj.mass = const.mRb;
            elseif strcmpi(atomType,'K40')
                obj.mass = const.mK;
            else
                error('Atom type not supported!');
            end
            obj.a = odt.scale(obj.wavelength,obj.atomType);
            if nargin < 3
                obj.useGravity = true;
            else
                obj.useGravity = logical(useGravity);
            end
        end
        
        function checkInputs(obj,P,wx,wy)
            if numel(wx)>1 && numel(P)>1
                error('Only one of wx and P can have more than one element!');
            end
            if nargin==4
                if numel(wy)>1 && (numel(wx)>1 || numel(P)>1)
                    error('Only one of wx, wy, and P can have more than one element!');
                end
                obj.wx = wx;
                obj.wy = wy;
            else
                obj.wx = wx;
                obj.wy = wx;
            end
            obj.P = P;
        end
        
        function E = U(obj,x,P,wx,wy)
            if nargin<=4
                wy = wx;
            end
            if obj.useGravity
                E = -2*obj.a*P./(pi*wx.*wy).*exp(-2*x.^2./wx.^2)-obj.mass*const.g*x;
            else
                E = -2*obj.a*P./(pi*wx.*wy).*exp(-2*x.^2./wx.^2);
            end
        end
        
        function [xMin,xMax] = minmax(obj,P,wx,wy)
            if nargin < 4
                wy = wx;
            end
            obj.checkInputs(P,wx,wy);
            
            if ~obj.useGravity
                if numel(obj.wx) > 1
                    xMax = Inf(size(obj.w));
                    xMin = zeros(size(obj.w));
                else
                    xMax = Inf(size(obj.P));
                    xMin = zeros(size(obj.P));
                end
            else
                wArg = -pi^2*obj.mass^2*const.g^2*obj.wx.^4.*obj.wy.^2./(16*obj.a.^2.*obj.P.^2);
                xMin = zeros(size(wArg));
                xMax = xMin;
                if numel(obj.wx) > 1
                    for nn=1:numel(obj.wx)
                        xMin(nn) = sqrt(-Lambert_W(wArg(nn),0)).*obj.wx(nn)/2;
                        xMax(nn) = sqrt(-Lambert_W(wArg(nn),-1))*obj.wx(nn)/2;
                    end
                else
                    for nn=1:numel(obj.P)
                        xMin(nn) = sqrt(-Lambert_W(wArg(nn),0)).*obj.wx/2;
                        xMax(nn) = sqrt(-Lambert_W(wArg(nn),-1))*obj.wx/2;
                    end
                end
            end
        end
        
        function f = freq(obj,P,wx,wy)
            if nargin<4
                wy = wx;
            end
            obj.checkInputs(P,wx,wy);
            xMin = obj.minmax(P,wx,wy);
            f = (1./(2*pi))*real(sqrt(8*obj.a*obj.P./(pi*obj.mass*obj.wx.^3*obj.wy).*(1-4*xMin.^2./obj.wx.^2).*exp(-2*xMin.^2./obj.wx.^2)));
        end
        
        function f = axialFreq(obj,P,wx,wy)
            if nargin<4
                wy = wx;
            end
            obj.checkInputs(P,wx,wy);
            f = (1./(2*pi)).*sqrt(4*obj.a.*obj.wavelength.^2.*obj.P./(pi^3*obj.mass.*obj.wx.^3*obj.wy^3));
        end
        
        function D = depth(obj,P,wx,wy)
            if nargin<4
                wy = wx;
            end
            obj.checkInputs(P,wx,wy);
            [xMin,xMax] = obj.minmax(P,wx,wy);
            D = zeros(size(xMin));
            for nn=1:numel(xMin)
                if numel(obj.wx) > 1
                    D(nn) = abs(obj.U(xMin(nn),obj.P,obj.wx(nn),obj.wy)-obj.U(xMax(nn),obj.P,obj.wx(nn),obj.wy));
                else
                    D(nn) = abs(obj.U(xMin(nn),obj.P(nn),obj.wx,obj.wy)-obj.U(xMax(nn),obj.P(nn),obj.wx,obj.wy));
                end
            end
        end
        
        function w0 = solveWaist(obj,P,f,wy)
            if ~obj.useGravity
                if nargin == 3
                    w0 = (8*obj.a*P./(pi*obj.mass*(2*pi*f).^2)).^(1/4);
                else
                    w0 = (8*obj.a*P./(pi*obj.mass*(2*pi*f).^2*wy)).^(1/3);
                end
            else
                options = optimset('display','off');
                if nargin == 3
                    w0 = fsolve(@(x) obj.freq(P,x) - f,60e-6,options);
                else
                    w0 = fsolve(@(x) obj.freq(P,x,wy) - f,60e-6,options);
                end
            end
        end
    end
    
    methods(Static)
        function a=scale(wavelength,atomType)
            if strcmpi(atomType,'Rb87')
                atom = odt.Rb87;
            elseif strcmpi(atomType,'K40')
                atom = odt.K40;
            else
                error('Atom type not supported!');
            end
            w = 2*pi*const.c/wavelength;
            a = 0;
            for nn=1:numel(atom)
                w0 = 2*pi*const.c/atom(nn).wavelength;
                a = a + (2*atom(nn).J+1)/2*pi*const.c^2/(2*w0.^3).*atom(nn).decay.*((w0+w).^(-1)+(w0-w).^(-1));
            end
        end
    end
    
end