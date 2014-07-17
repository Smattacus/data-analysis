function [ph] = genPhase(s1, center, delta, dt)
%Function which generates the angle array of the modulation frequency
%from an array of PMT readings.
%[phase]] = genPhase(time_series, chop, chop_delta, dt);
%INPUTS:
%time_series    = Data in time.
%chop           = Frequency of the chop.
%chop_delta     = 1/2 the Width of the window to cut off with.
%dt             = Time between points on time_series.
[f, g] = spec(s1, dt);
ig = find(abs(abs(f) - center) > delta);
g(ig) = 0;
[t, x] = ispec(g, f);
ph = unwrap(angle(hilbert(x)));
