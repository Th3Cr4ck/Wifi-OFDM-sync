clear
close all
clc

addpath('..');

pilotValues = [1 1 1 1];
channelGainAtPilots = exp(1j*pi/4) * ones(1,4);

residualCfoPhase = estimateResidualCfoPhase(pilotValues, channelGainAtPilots);

expected = -pi/4;

fprintf("Estimado: %.4f rad\n", residualCfoPhase);
fprintf("Esperado: %.4f rad\n", expected);

assert(abs(residualCfoPhase - expected) < 1e-10, ...
    "ERROR: resultado incorrecto.");

fprintf("TEST OK\n");