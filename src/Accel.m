%--------------------------------------------------------------------------
%
% Accel.m
%
% Purpose:
%   Computes the acceleration of an Earth orbiting satellite due to 
%    - the Earth's harmonic gravity field, 
%    - the gravitational perturbations of the Sun and Moon
%    - the solar radiation pressure and
%    - the atmospheric drag
%
% Inputs:
%   Mjd_TT      Terrestrial Time (Modified Julian Date)
%   Y           Satellite state vector in the ICRF/EME2000 system
%
% Output:
%   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function [dY] = Accel(x, Y)

global AuxParam eopdata

SAT_Const

[UT1_UTC, TAI_UTC, x_pole, y_pole] = IERS(eopdata, AuxParam.Mjd_TT + x/86400);
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);

Mjd_UT1 = AuxParam.Mjd_TT + x/86400 + (UT1_UTC-TT_UTC)/86400.0;

P = PrecMatrix(MJD_J2000,AuxParam.Mjd_TT + x/86400);
N = NutMatrix(AuxParam.Mjd_TT + x/86400);
T = N * P;
E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

% Acceleration due to harmonic gravity field
a = AccelHarmonic(Y(1:3), E, AuxParam.n, AuxParam.m);

dY = [Y(4:6);a];

