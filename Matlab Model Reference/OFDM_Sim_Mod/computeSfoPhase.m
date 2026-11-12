function sfoPhase = computeSfoPhase(OFDMSymbolNumber, Nfft, Ng,...
                                    prevEstSfo, activeIndex)
% COMPUTESFOPHASE Computes the phase correction due to SFO for the
% l'th OFDM Symbol.
%
% Input Arguments:
%  - OFDMSymbolNumber: OFDM symbol index l
%  - Nfft: FFT size
%  - Ng: Guard samples length
%  - prevEstSfo: Previous OFDM symbol SFO estimation
%  - activeIndex: Active subcarrier indexes
%
% Output Arguments:
%  - sfoPhase: Phase correction for each active subcarrier.

    % Active subcarrier numbers
    % k = [-26:-1 1:26];

    sfoPhase = 2*pi*OFDMSymbolNumber * (Nfft+Ng)/Nfft * prevEstSfo .* activeIndex;
end