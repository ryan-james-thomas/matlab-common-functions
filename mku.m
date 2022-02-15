classdef mku < handle
    %MKU Class definition for handling the control of the MKU PLL
    %synthesizer via an ESP8266 WiFi micro-controller
    %
    properties(SetAccess = protected)
        conn    %TCPIP object handling the connection
        timeout %Timeout for reading data
    end
    
    properties
        host    %IP address of ESP8266 WiFi server
        port    %Port to connect to
        f1      %F1 frequency in Hz
        f2      %F2 frequency in Hz
    end
    
    properties(Constant)
        DEFAULT_HOST = '192.168.1.33';
        DEFAULT_PORT = 6666;
        DEFAULT_TIMEOUT = 1;
    end
    
    methods
        function self = mku(host,port)
            %MKU Creates the MKU object.
            %
            %   SELF = MKU() Creates MKU with default host name and port
            %   number
            %
            %   SELF = MKU(HOST) Creates MKU with given host name and
            %   default port number
            %
            %   SELF = MKU(HOST,PORT) Creates MKU with given host name and
            %   port number
            %
            self.host = self.DEFAULT_HOST;
            self.port = self.DEFAULT_PORT;
            self.timeout = self.DEFAULT_TIMEOUT;
            if nargin == 1
                self.host = host;
            elseif nargin == 2
                self.host = host;
                self.port = port;
            end
        end
        
        function self = setTimeout(self,timeout)
            %SETTIMEOUT Sets the timeout for reading data from the ESP8266
            %device
            %
            %   SELF = MKU.SETTIMEOUT(TIMEOUT) Sets the timeout to TIMEOUT,
            %   in seconds
            %
            self.timeout = timeout;
        end
        
        function self = open(self)
            %OPEN Opens a connection to the host
            r = instrfindall('RemoteHost',self.host,'RemotePort',self.port);
            if isempty(r)
                %
                % If no connections exist, create that connection
                %
                self.conn = tcpip(self.host,self.port);
                self.conn.InputBufferSize = 2^24;
                self.conn.OutputBufferSize = 2^24;
                fopen(self.conn);
            elseif strcmpi(r.Status,'closed')
                %
                % If a connection exists but it is closed, set the buffer
                % size correctly and then open it
                %
                self.conn = r;
                self.conn.InputBufferSize = 2^24;
                self.conn.OutputBufferSize = 2^24;
                fopen(self.conn);
            else
                %
                % Otherwise set the client parameter to that connection
                %
                self.conn = r;
            end
        end
        
        function self = close(self)
            %CLOSE Closes the connection to the host
            if ~isempty(self.conn) && isvalid(self.conn) && strcmpi(self.conn.Status,'open')
                fclose(self.conn);
            end
            delete(self.conn);
            self.conn = [];
        end
        
        function delete(self)
            %DELETE Deletes the current object. Closes the connection first
            try
                self.close;
            catch
                disp('Error deleting client');
            end
        end
        
        function r = cmd(self,cmd)
            %CMD Writes a command to the device and returns the response
            %
            %   R = MKU.CMD(CMD) Writes the command CMD (a character
            %   vector) and returns the response R (a character vector)
            %
            self.open;
            fprintf(self.conn,cmd);
            pause(10e-3);
            jj = 1;
            while ~self.conn.BytesAvailable
                pause(10e-3);
                if jj > floor(self.timeout/10e-3)
                    error('Timeout reading data');
                end
                jj = jj + 1;
            end
            pause(100e-3);
            r = char(fread(self.conn,self.conn.BytesAvailable))';
            self.close;
        end
        
        function self = set(self,f1,f2)
            if nargin == 2
                if numel(f1) == 1
                    self.f1 = f1;
                elseif numel(f1) == 2
                    self.f1 = f1(1);
                    self.f2 = f1(2);
                end
            elseif nargin == 3
                self.f1 = f1;
                self.f2 = f2;
            end
        end
        
        function r = status(self)
            r = self.cmd('sa');
        end
        
        function self = writeList(self,varargin)
            self.set(varargin{:});
            self.check(self.f1);
            self.check(self.f2);
            
            r = self.cmd(self.parseFreq(self.f1,1));
            if any(r == 'N')
                error('Error writing frequency 1: %s',strrep(r,newline,''));
            end
            r = self.cmd(self.parseFreq(self.f2,2));
            if any(r == 'N')
                error('Error writing frequency 2: %s',strrep(r,newline,''));
            end
        end
        
        
    end
    
    methods(Static)
        function s = parseFreq(f,ch)
            if nargin < 2
                ch = 1;
            elseif ch ~= 1 && ch ~= 2
                error('Channel can only be 1 or 2')
            end
                
            fG = floor(f/1e9);
            fM = floor((f - fG*1e9)/1e6);
            fK = floor((f - fG*1e9 - fM*1e6)/1e3);
            fH = floor((f - fG*1e9 - fM*1e6 - fK*1e3));
            s = sprintf('%03dGF%d%03dMF%d%03dkF%d%03dHF%d',fG,ch,fM,ch,fK,ch,fH,ch);
        end
        
        function check(f)
            if (f < 54e6 || f > 6850e6) && (f < 8e9 || f > 13e9)
                error('Frequency must be between [54,6850] MHz (Aux) or [8,13] GHz (Main)');
            end
        end
    end
end

