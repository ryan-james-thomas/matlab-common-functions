classdef RS_Synthesiser < handle
    properties(SetAccess = protected)
        dev
        visa_label
    end
    
    properties
        mode
        freq
        pow
        list
        state
    end
    
    properties(Constant)
        DEFAULT_VISA_LABEL = 'TCPIP::192.168.1.11::INSTR';
    end
    
    methods
        function self = RS_Synthesiser(visa_label)
            if nargin < 1
                self.visa_label = RS_Synthesiser.DEFAULT_VISA_LABEL;
            else
                self.visa_label = visa_label;
            end
            
            self.open;
        end
        
        function self = open(self)
            if isempty(self.dev) || ~isvalid(self.dev)
                self.dev = visa('rs',self.visa_label);
            end
            
            if strcmpi(self.dev.Status,'closed')
                fopen(self.dev);
            end
        end
        
        function self = close(self)
            if ~isempty(self.dev) && isvalid(self.dev)
                if strcmpi(self.dev.Status,'open')
                    fclose(self.dev);
                end
                delete(self.dev);
                self.dev = [];
            end
        end
        
        function self = write(self,cmd,varargin)
            s = sprintf('%s\n',sprintf(cmd,varargin{:}));
            self.open;
            fprintf(self.dev,s);
        end
        
        function s = read(self)
            s = fgetl(self.dev);
        end
        
        function self = writeState(self,s)
            if ischar(s)
                if strcmpi(s,'on')
                    self.state = 'on';
                elseif strcmpi(s,'off')
                    self.state = 'off';
                end
            else
                if s
                    self.state = 'on';
                else
                    self.state = 'off';
                end
            end 
            self.write('output:state %s',self.state);
        end
        
        function self = writeCW(self,f,p)
            if nargin == 2
                self.freq = f;
            end
            if nargin == 3
                self.pow = p;
            end
            self.mode = 'cw';
            
            self.write('source:freq:freq %.6f',self.freq);
            self.write('source:pow:pow %.3f',self.pow);
            self.write('source:freq:mode %s',self.mode);
            
        end
        
        function self = writeList(self,f,p)
            if nargin >= 2
                self.list.freq = f;
            end
            if nargin >= 3
                self.list.pow = p;
            end
            self.mode = 'list';
            
            self.write('source:list:set "New_List"');
            self.write('source:list:dwell 1ms');
            self.write('source:list:trig:source external');
            self.write('source:list:mode step');
            f = num2cell(self.list.freq);
            s = sprintf('%.6f MHz,',f{:});
            s = s(1:end-1);
            self.write('source:list:freq %s',s);
            p = num2cell(self.list.pow);
            s = sprintf('%.6f dBm,',p{:});
            s = s(1:end-1);
            self.write('source:list:pow %s',s);
            self.write('source:freq:mode %s',self.mode);
        end
        
    end
    
end