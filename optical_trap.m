classdef optical_trap < handle
    properties
        atom_type
        mass
        lasers
        ext_force
    end
    
    properties(Constant,Hidden = true)
        Rb87 = [struct('wavelength',780e-9,'decay',2*pi*6.065e6,'J',3/2);
                struct('wavelength',795e-9,'decay',2*pi*5.74e6,'J',1/2)];
            
        K40 = [struct('wavelength',767e-9,'decay',2*pi*6.035e6,'J',3/2);
               struct('wavelength',770e-9,'decay',2*pi*5.956e6,'J',1/2)];

    end
    
    methods
        function self = optical_trap(atom_type,lasers)
            if strcmpi(atom_type,'Rb87')
                self.atom_type = 'Rb87';
                self.mass = const.mRb;
            elseif strcmpi(atom_type,'K40')
                self.atom_type = 'K40';
                self.mass = const.mK;
            else
                error('Only Rb87 and K40 are currently supported!');
            end
            
            if nargin > 1 && isa(lasers,'gaussian_beam')
                self.lasers = lasers;
            end
                
        end
        
        function U = potential(self,x,y,z)
            U = 0;
            for nn = 1:numel(self.lasers)
                U = U - self.scale(self.lasers(nn).wavelength,self.atom_type)*self.lasers(nn).intensity(x,y,z);
            end
            
        end
        
        function [Fx,Fy,Fz] = force(self,x,y,z)
            [Fx,Fy,Fz] = deal(0);
            for nn = 1:numel(self.lasers)
                [tmpx,tmpy,tmpz] = self.lasers(nn).force(x,y,z);
                Fx = Fx + self.scale(self.lasers(nn).wavelength,self.atom_type)*tmpx;
                Fy = Fy + self.scale(self.lasers(nn).wavelength,self.atom_type)*tmpy;
                Fz = Fz + self.scale(self.lasers(nn).wavelength,self.atom_type)*tmpz;
            end
            
            if ~isempty(self.ext_force)
                Fx = Fx + self.ext_force.Fx(x,y,z);
                Fy = Fy + self.ext_force.Fy(x,y,z);
                Fz = Fz + self.ext_force.Fz(x,y,z);
            end
        end
        
        function [f,V,K,r0] = freq(self,x,y,z,dr,opt)
            K = zeros(3,3);
            if nargin < 6
                opt = 1;
            end
            if opt
                r0 = self.find_zero(x,y,z);
                x = r0(1);y = r0(2);z = r0(3);
            end
            args = {x + dr*[0,1],y,z;
                    x,y + dr*[0,1],z;
                    x,y,z + dr*[0,1]};
            for row = 1:3
                [Fx,Fy,Fz] = self.force(args{row,:});    
                
                K(row,1) = (Fx(2) - Fx(1))/dr;
                K(row,2) = (Fy(2) - Fy(1))/dr;
                K(row,3) = (Fz(2) - Fz(1))/dr;
            end
            
            [V,D] = eig(K);
            f = sqrt(diag(D)/const.mRb)/(2*pi);
        end
        
        function r0 = find_zero(self,x,y,z)
            function F = int_force(r)
                [Fx,Fy,Fz] = self.force(r(1),r(2),r(3));
                F = [Fx;Fy;Fz];
            end
            
            r0 = fsolve(@(x) int_force(x)/self.mass,[x;y;z],optimset('display','off'));
                
        end
        
        function isosurface(self,X,Y,Z,isovalue)
            U = self.potential(X,Y,Z)/const.kb*1e6;
            fv = isosurface(X,Y,Z,U,isovalue);
            p = patch(fv);
            isonormals(X,Y,Z,U,p);
            p.FaceColor = 'red';
            p.EdgeColor = 'none';
            daspect([1,1,1]);
            view(3);
            axis tight;
            camlight;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
        end
    end
    
    methods(Static)
        function a = scale(wavelength,atom_type)
            if strcmpi(atom_type,'Rb87')
                atom = optical_trap.Rb87;
            elseif strcmpi(atom_type,'K40')
                atom = optical_trap.K40;
            else
                error('Atom type not supported!');
            end
            w = 2*pi*const.c/wavelength;
            a = 0;
            for nn = 1:numel(atom)
                w0 = 2*pi*const.c/atom(nn).wavelength;
                a = a + (2*atom(nn).J+1)/2*pi*const.c^2/(2*w0.^3).*atom(nn).decay.*((w0+w).^(-1)+(w0-w).^(-1));
            end
        end
    end
    
    
end