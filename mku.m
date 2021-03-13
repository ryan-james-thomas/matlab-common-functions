classdef mku < handle
%MKU Defines a class that interfaces with the Kuhne MKU PLL synthesizer via
%a NodeMCU module
%
% M = MKU(F1,F2) with optional arguments F1 and F2 sets the internal
% properties F1 and F2 to their respective values.  These are the output
% frequencies of the MKU.  To write to the device, use the WRITEF1 and
% WRITEF2 methods
    properties
        f1      %Frequency 1, output when SW1 is open or high
        f2      %Frequency 2, output when SW1 is closed or low
        msg     %Message returned from NodeMCU interface
    end
    
    properties(Constant)
        uri = 'http://172.22.251.245/mku/'  %URI for the MKU
    end
    
    methods
        function obj = mku(f1,f2)
            if nargin == 1
                mku.setf1(f1);
            elseif nargin == 2
                mku.setf1(f1);
                mku.setf2(f2);
            end
        end
        
        %SETF1 Sets the internal frequency F1 to the desired value
        %
        %OBJ = SETF1(F) Sets the internal frequency F1 to f
        function obj = setf1(obj,f)
            obj.check(f);
            obj.f1 = f;
        end
        
        %SETF2 Sets the internal frequency F2 to the desired value
        %
        %OBJ = SETF2(F) Sets the internal frequency F2 to f
        function obj = setf2(obj,f)
            obj.check(f);
            obj.f2 = f;
        end
        
        %WRITEF1 Writes internal frequency F1 to the MKU
        %
        %OBJ = WRITEF1() Writes internal frequency F1 to the MKU
        function obj = writef1(obj)
            obj.check(obj.f1);
            q = obj.format(obj.f1);
            r = webread([obj.uri,'F1'],q{:});
            obj.msg = jsondecode(r);
        end
        
        %WRITEF2 Writes internal frequency F2 to the MKU
        %
        %OBJ = WRITEF2() Writes internal frequency F2 to the MKU
        function obj = writef2(obj)
            obj.check(obj.f2);
            q = obj.format(obj.f2);
            r = webread([obj.uri,'F2'],q{:});
            obj.msg = jsondecode(r);
        end
        
    end
    
    methods(Static)
        %CHECK Checks that the frequency is in the right range
        %
        %CHECK(F) Checks that F is within the MKU range of [54,6850] MHz
        %or [8,13] GHz
        function check(f)
            if (f<54 || f>6850) && (f<8e3 || f>13e3)
                error('Frequency must be between [54,6850] MHz (Aux) or [8,13] GHz (Main)');
            end
        end
        
        %FORMAT Converts a frequency to a cell array of arguments
        %
        %S = FORMAT(F) Converts frequency F to a cell array of values that
        %can be interpreted by the NodeMCU firmware and is easily written
        %to the device using the WEBREAD MATLAB function
        function s = format(f)
            G = floor(f/1e3);
            M = floor((f-G*1e3));
            k = floor((f-G*1e3-M)*1e3);
            H = floor((f-G*1e3-M-k*1e-3)*1e6);
            s = {'G',sprintf('%03d',G),...
                                'M',sprintf('%03d',M),...
                                'k',sprintf('%03d',k),...
                                'H',sprintf('%03d',H)};
        end
        
    end
    
end