classdef optical_trap < handle
    properties
        atom_type
        mass
        lasers
        external_acceleration
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