
function [ZoneBeta ZoneLambda time Mag Location LogPSA] = MonteCarloPSHA(Mw,Years,MinMw,MaxMw,ZoneBoundary,SeismicBoundary,NumZones,Site,GMPEConsidered,GMPEWeights,SimPeriod,NumRuns)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[ZoneBeta ZoneNo ZoneLambda ] = RecurrenceRate(Mw,Years,3.0,7)
[time Mag Location ] = SyntheticCat( ZoneLambda,ZoneBeta,MinMw,MaxMw,ZoneBoundary,SeismicBoundary,NumZones,SimPeriod,NumRuns)
[LogPSA ] = GMPEEstimates(Mag,Location,Site,GMPEConsidered,GMPEWeights)


end

