function out=bin2hex(str)
    out = dec2hex(bin2dec(str),length(str)/4);
end