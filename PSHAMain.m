clear;
clc;
load GMPE_data.mat
load SeismicBoundary.mat
load Mw.mat
load Years.mat
load ZoneBoundary.mat
[ZoneBeta ZoneLambda time Mag Location LogPSA] = MonteCarloPSHA(Mw,Years,4.0,7.5,ZoneBoundary,SeismicBoundary,1,[0.2,5.61],{'AtkinsonHardrock';'AtkinsonBC';'PezeshkHybridEM';'SilvaDCM';'SilvaDCMSaturation';'SilvaScmConStressSatu';'SilvaScmVarStress';'SilvaScmConstStress'},[0.1;0.2;0.1;0.1;0.1;0.2;0.1;0.1],100,1);