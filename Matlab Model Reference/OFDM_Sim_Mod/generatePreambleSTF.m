function preamShort = generatePreambleSTF()
    k    = [-24 -20 -16 -12 -8 -4 4 8 12 16 20 24];
    vals = (1+1j) * [1 -1 1 -1 -1 1 -1 -1 1 1 1 1];
    
    S = zeros(1,64);
    S(mod(k,64)+1) = sqrt(13/6) * vals;   % mapeo a bins FFT (k negativo -> 64+k)
    
    s64   = ifft(S);                       % 64 muestras, periodo de 16
    preamShort = repmat(s64(1:16), 1, 10);      % 160 muestras = 8 µs
end