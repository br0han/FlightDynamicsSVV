clear; close all; clc;
run('flightdataAnalysis.m');
close all;

[P,Periods] = Period(D_Dutch_yawrate,1,1e-3) %fill in what (t,y) array to calculate period add expected distance between peaks


%Period Calculator (uses peaks and the distance between these as period)
function [Q,Periods] = Period(Phugoid_V,dist,Thres)
    t = Phugoid_V(:,1); y = Phugoid_V(:,2);
    y = max(max(y)-y,y-min(y));
    [PV_y,PV_t] = findpeaks(y,t,'minpeakdistance',dist,'Threshold',Thres);
    peaks = [];
    for i = 1:length(PV_t)
        t_i = PV_t(i);
        num = find(Phugoid_V==t_i);
        x = Phugoid_V(num,2);
        peaks = [peaks ; x];
    end
    PV_t = [PV_t, peaks];
    final_V = mean(PV_t(:,2));
    Periods = [];
    for i = 1:length(PV_t)-1
        Per  = 2*(PV_t(i+1,1)-PV_t(i,1));
        Periods = [Periods ; Per];
        if i ~= length(PV_t)-1
            Per = (PV_t(6,1)-PV_t(i,1))/((-i+6)/2);
            Periods = [Periods ; Per];
        end
    end
    Q = mean(Periods);
end

    