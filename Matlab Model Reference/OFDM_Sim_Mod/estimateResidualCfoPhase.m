function residualCfoPhase = estimateResidualCfoPhase(pilotValues, channelGainAtPilots)
% ESTIMATERESIDUALCFOPHASE Computes the estimated residual CFO based on the
% channel gains at the pilot subcarriers.
% 
% Input Arguments:
%   - pilotValues: Vector of the values of the pilots.
%   - channelGainAtPilots: Vector of the estimated channel gains at pilot subcarriers.
%
% Output Arguments:
%   - residualCfo: Estimated residual CFO.


correlation = sum(pilotValues .* conj(channelGainAtPilots));
residualCfoPhase = angle(correlation);

end