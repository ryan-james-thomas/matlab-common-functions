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

            self.ext_force.Fx = @(x,y,z) 0;
            self.ext_force.Fy = @(x,y,z) 0;
            self.ext_force.Fz = @(x,y,z) 0;
                
        end

        function self = set_powers(self,P)
            if numel(P) == 1
                P = P*ones(numel(self.lasers),1);
            end
            for nn = 1:numel(self.lasers)
                self.lasers(nn).power = P(nn);
            end
        end
        
        function U = potential(self,x,y,z)
            U = 0;
            for nn = 1:numel(self.lasers)
                U = U - self.scale(self.lasers(nn).wavelength,self.atom_type)*self.lasers(nn).intensity(x,y,z);
            end
        end
        
        function [Fx,Fy,Fz] = force(self,x,y,z,opt)
            [Fx,Fy,Fz] = deal(0);
            for nn = 1:numel(self.lasers)
                [tmpx,tmpy,tmpz] = self.lasers(nn).force(x,y,z);
                Fx = Fx - self.scale(self.lasers(nn).wavelength,self.atom_type)*tmpx;
                Fy = Fy - self.scale(self.lasers(nn).wavelength,self.atom_type)*tmpy;
                Fz = Fz - self.scale(self.lasers(nn).wavelength,self.atom_type)*tmpz;
            end
            
            if ~isempty(self.ext_force)
                Fx = Fx + self.ext_force.Fx(x,y,z);
                Fy = Fy + self.ext_force.Fy(x,y,z);
                Fz = Fz + self.ext_force.Fz(x,y,z);
            end

            if nargin > 4
                switch lower(opt)
                    case 'y'
                        Fx = Fy;
                    case 'z'
                        Fx = Fz;
                end
            end
        end
        
        function [f,V,K,r0] = freq(self,x,y,z,dr,opt)
            K = zeros(3,3);
            if nargin < 6
                opt = 1;
            end
            if opt
                r0 = self.find_zero(x,y,z);
                if isempty(r0)
                    f = zeros(3,1);
                    V = K;
                    return
                end
                x = r0(1);y = r0(2);z = r0(3);
            end
            args = {x + dr*[0,1],y,z;
                    x,y + dr*[0,1],z;
                    x,y,z + dr*[0,1]};
            for row = 1:3
                [Fx,Fy,Fz] = self.force(args{row,:});    
                
                K(row,1) = -(Fx(2) - Fx(1))/dr;
                K(row,2) = -(Fy(2) - Fy(1))/dr;
                K(row,3) = -(Fz(2) - Fz(1))/dr;
            end
            
            [V,D] = eig(K);
            ftmp = sqrt(diag(D)/const.mRb)/(2*pi);
            f = zeros(size(ftmp));
            idx = zeros(size(ftmp));
            for nn = 1:numel(f)
                [~,idx(nn)] = max(abs(V(:,nn).^2));
                f(idx(nn)) = ftmp(nn);
            end
            V = V(:,idx);

        end
        
        function r0 = find_zero(self,x,y,z,force_output)
            if nargin < 5
                force_output = false;
            end

            function F = int_force(r)
                [Fx,Fy,Fz] = self.force(r(1),r(2),r(3));
                F = [Fx;Fy;Fz];
            end
            
            r0 = fsolve(@(x) int_force(x)/self.mass,[x;y;z],optimset('display','off'));
            Fx = self.force(r0(1) + [0,1e-6],r0(2),r0(3),'x');
            Fy = self.force(r0(1),r0(2) + [0,1e-6],r0(3),'y');
            Fz = self.force(r0(1),r0(2),r0(3) + [0,1e-6],'z');
            if (diff(Fx) > 0 || diff(Fy) > 0 || diff(Fz) > 0) && ~force_output
                r0 = [];
            end

        end

        function D = depth(self,x,y,z,direction)
            r0 = self.find_zero(x,y,z);
            if isempty(r0)
                D = 0;
                return
            end
            r1 = r0;
            dr = 10e-6;
            switch lower(direction)
                case 'x'
                    F = self.force(r0(1) + [0,1e-7],r0(2),r0(3));
                    if diff(F) > 0
                        xs = r0(1) + 1e-7;
                        result = search_for_zero(@(x) self.force(x,r0(2),r0(3),'x'),xs,5e-3,dr);
                    else
                        xs = r0(1) - 1e-7;
                        result = search_for_zero(@(x) self.force(x,r0(2),r0(3),'x'),xs,-5e-3,dr);
                    end
                    r1(1) = result;
                    xx = linspace(r0(1),r1(1),1e2);
                    D = trapz(xx,self.force(xx,r0(2),r0(3),'x'));

                case 'y'
                    F = self.force(r0(1),r0(2) + [0,1e-7],r0(3),'y');
                    if diff(F) > 0
                        xs = r0(2) + 1e-7;
                        result = search_for_zero(@(x) self.force(r0(1),x,r0(3),'y'),xs,5e-3,dr);
                    else
                        xs = r0(2) - 1e-7;
                        result = search_for_zero(@(x) self.force(r0(1),x,r0(3),'y'),xs,-5e-3,dr);
                    end
                    r1(2) = result;
                    xx = linspace(r0(2),r1(2),1e2);
                    D = trapz(xx,self.force(r0(1),xx,r0(3),'y'));

                case 'z'
                    F = self.force(r0(1),r0(2),r0(3) + [0,1e-7],'z');
                    if diff(F) > 0
                        xs = r0(3) + 1e-7;
                        result = search_for_zero(@(x) self.force(r0(1),r0(2),x,'z'),xs,5e-3,dr);
                    else
                        xs = r0(3) - 1e-7;
                        result = search_for_zero(@(x) self.force(r0(1),r0(2),x,'z'),xs,-5e-3,dr);
                    end
                    r1(3) = result;
                    xx = linspace(r0(3),r1(3),1e2);
                    D = trapz(xx,self.force(r0(1),r0(2),xx,'z'));

            end

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

function r = search_for_zero(func,xs,xmax,dx)
    dx = dx*sign(xmax - xs);
    x0 = xs;
    while sign(func(x0)) == sign(func(xs))
        xs = x0;
        x0 = x0 + dx;
        if abs(x0) > abs(xmax)
            r = xmax;
            return;
        end
    end

    r = fsolve(func,0.5*(x0 + xs),optimset('display','off'));

end