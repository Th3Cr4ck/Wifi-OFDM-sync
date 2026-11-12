function [fineCfo, estOffset, W] = estimateFineOffset(pilotValues, ...
    prevPilotValues, Tsym, fc, initialOffset, W)
% ESTIMATEFINEOFFSET Estimates the fine frequency offset using pilot values.
%
% Input Arguments:
%   - pilotValues:      Pilot values of the current OFDM symbol.
%   - prevPilotValues:  Pilot values of the previous OFDM symbol.
%   - Tsym:             OFDM symbol duration.
%   - fc:               Carrier frequency.
%   - initialOffset:    Initial frequency-offset estimate (by S&C, eps_0).
%
% Output Arguments:
%   - fineCfo:           Estimated residual CFO (eps_r).
%   - estOffset:         Updated frequency-offset estimate (eps=eps_0+eps_r).

% W_l = sum_i P_{l,i} * conj(P_{l-1,i}), acumulada entre simbolos.
% fineCfo estima eps_r; estOffset = eps_0 + eps_r estima eps.
% initialOffset es el eps_0 fijo (LT), no la estimacion anterior.

% Calculate W_l (autocorrelation)
Wl = sum(pilotValues .* conj(prevPilotValues));

% Continuously accumulate W_l
W = W + Wl;

% Estimate residual CFO
fineCfo = angle(W) / (2*pi*Tsym*fc);

% Update total frequency-offset estimate
estOffset = initialOffset + fineCfo;

end