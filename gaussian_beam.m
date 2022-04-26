classdef gaussian_beam < handle
    properties
        wavelength
        power
        waists
        center
        rotation_angles     
    end
    
    properties(SetAccess = protected)
        rotation_matrix
    end
    
    methods
        function self = gaussian_beam(wavelength,P,w,r0,rotation_angles)
            self.wavelength = wavelength;
            self.power = P;
            if numel(w) == 1
                self.waists = w*[1,1];
            elseif numel(w) == 2
                self.waists = w(:)';
            else
                error('Only two waist values can be given!');
            end
            %
            % Set center positions
            %
            if nargin > 3
                if numel(r0) == 1 
                    self.center = r0*[1,1,1];
                elseif numel(r0) == 3
                    self.center = r0(:)';
                else
                    error('Only three center positions (x,y,z) can be given!');
                end
            else
                self.center = [0,0,0];
            end
            %
            % Set rotation angles and matrix
            %
            if nargin > 4
                self.set_rotations(rotation_angles)
            else
                self.set_rotations(0);
            end
            
        end
        
        function self = set_rotations(self,rotation_angles)
            self.rotation_angles = rotation_angles;
            if numel(self.rotation_angles) == 1
                th = self.rotation_angles(1);
                self.rotation_matrix = [cosd(th), 0, sind(th);
                                        0,  1,  0;
                                        -sind(th),  0,  cosd(th)];
            elseif numel(self.rotation_angles) == 2
                th = self.rotation_angles(1);
                R1 = [cosd(th),     sind(th),   0;
                      -sind(th),    cosd(th),   0;
                      0,            0,          1];
                th = self.rotation_angles(2);
                R2 = [cosd(th), 0, sind(th);
                      0,  1,  0;
                      -sind(th),  0,  cosd(th)];
                self.rotation_matrix = R2*R1;
            elseif numel(self.rotation_angles) == 3
                th = self.rotation_angles(1);
                R1 = [cosd(th),     sind(th),   0;
                      -sind(th),    cosd(th),   0;
                      0,            0,          1];
                th = self.rotation_angles(2);
                R2 = [cosd(th), 0, sind(th);
                      0,  1,  0;
                      -sind(th),  0,  cosd(th)];
                th = self.rotation_angles(3);
                R3 = [1,    0,      0;
                      0,    cosd(th),   sind(th);
                      0,    -sind(th),  cosd(th)];
                self.rotation_matrix = R3*R2*R1;
            end
        end
        
        function I = intensity(self,x,y,z)
            x_lab = x - self.center(1);
            y_lab = y - self.center(2);
            z_lab = z - self.center(3);
            
            R = self.rotation_matrix;
            %R takes body coordinates to lab coordinates
            x_body = R(1,1)*x_lab + R(2,1)*y_lab + R(3,1)*z_lab;
            y_body = R(1,2)*x_lab + R(2,2)*y_lab + R(3,2)*z_lab;
            z_body = R(1,3)*x_lab + R(2,3)*y_lab + R(3,3)*z_lab;

            zR = pi*self.waists.^2/self.wavelength;
            w_x = self.waists(1)*sqrt(1 + z_body.^2/zR(1)^2);
            w_y = self.waists(2)*sqrt(1 + z_body.^2/zR(2)^2);
            
            I = 2*self.power./(pi*w_x.*w_y).*exp(-2*x_body.^2./w_x.^2 - 2*y_body.^2./w_y.^2);
        end
        
        function [Fx_lab,Fy_lab,Fz_lab] = force(self,x,y,z)
            x_lab = x - self.center(1);
            y_lab = y - self.center(2);
            z_lab = z - self.center(3);
            
            R = self.rotation_matrix;
            %R takes body coordinates to lab coordinates
            x_body = R(1,1)*x_lab + R(2,1)*y_lab + R(3,1)*z_lab;
            y_body = R(1,2)*x_lab + R(2,2)*y_lab + R(3,2)*z_lab;
            z_body = R(1,3)*x_lab + R(2,3)*y_lab + R(3,3)*z_lab;

            zR = pi*self.waists.^2/self.wavelength;
            w_x = self.waists(1)*sqrt(1 + z_body.^2/zR(1)^2);
            w_y = self.waists(2)*sqrt(1 + z_body.^2/zR(2)^2);
            
            I = 2*self.power./(pi*w_x.*w_y).*exp(-2*x_body.^2./w_x.^2 - 2*y_body.^2./w_y.^2);
            Fx_body = 4*x_body./w_x.^2.*I;
            Fy_body = 4*y_body./w_x.^2.*I;
            Fz_body = -I./w_x.*(1 - 4*x_body.^2./w_x.^2).*self.waists(1)^2.*z_body./(zR(1)^2*w_x)...
                      -I./w_y.*(1 - 4*y_body.^2./w_y.^2).*self.waists(2)^2.*z_body./(zR(2)^2*w_y);
            
            Fx_lab = R(1,1)*Fx_body + R(1,2)*Fy_body + R(1,3)*Fz_body;
            Fy_lab = R(2,1)*Fx_body + R(2,2)*Fy_body + R(2,3)*Fz_body;
            Fz_lab = R(3,1)*Fx_body + R(3,2)*Fy_body + R(3,3)*Fz_body;
        end
        
    end
    
end