function out=hex2bin(str)
    out = dec2bin(hex2dec(str),length(str)*4);
end