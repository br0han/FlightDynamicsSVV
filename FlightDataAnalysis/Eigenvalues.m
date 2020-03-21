clear; close all; clc;
run('flightdataAnalysis.m');
close all;

dist = 1; 
Thres = 1e-3; 
%Period Calculator (uses peaks and the distance between these as period)

lst_dst = [10,10,1,1,1,1]; %minimum distance between peaks
lst_tres = [1e-6,1e-6,1e-4,1e-3,5e-4,4e-3]; %minimum hight change between peaks

%lst_data = [Phugoid_V,Phugoid_th,Dutch_yawrate,Dutch_rollrate,D_Dutch_yawrate,D_Dutch_rollrate];
s = 6;%choose which data set

if s == 1
    data = Phugoid_V;
end
if s == 2
    data = Phugoid_th;
end
if s == 3
    data = Dutch_yawrate;
end
if s == 4
    data = Dutch_rollrate;
end
if s == 5
    data = D_Dutch_yawrate;
end
if s == 6
    data = D_Dutch_rollrate;
end

Thres = lst_tres(s);
dist = lst_dst(s);



t = data(:,1); y = data(:,2);
y = max(max(y)-y,y-min(y));
[PV_y,PV_t] = findpeaks(y,t,'minpeakdistance',dist,'Threshold',Thres);
peaks = [];
for i = 1:length(PV_t)
    t_i = PV_t(i);
    num = find(data==t_i);
    x = data(num,2);
    peaks = [peaks ; x];
end
PV_t = [PV_t, peaks];
final = mean(PV_t(:,2));
Periods = [];
for i = 1:length(PV_t)-1
    Per  = 2*(PV_t(i+1,1)-PV_t(i,1));
    Periods = [Periods ; Per];
    if i ~= length(PV_t)-1
        Per = (PV_t(length(PV_t),1)-PV_t(i,1))/((-i+length(PV_t))/2);
        Periods = [Periods ; Per];
    end
end
Q = mean(Periods);


%Amplitudes
if s == 3 || s == 4 || s == 5 || s == 6
    final = data(1,2);
end

dev_avg = abs(data(:,2)-final);


figure();
title("peaks")
plot(data(:,1),data(:,2));
hold on
scatter(PV_t(:,1),PV_t(:,2),[],'r');
hold off

figure();
title("amplitudes")
plot(data(:,1),dev_avg);
hold on
scatter(PV_t(:,1),abs(PV_t(:,2)-final),[],'r');
hold off
    