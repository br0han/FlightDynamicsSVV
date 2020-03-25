clear; close all; clc;
run('flightdataAnalysis.m');
close all;

%Period Calculator (uses peaks and the distance between these as period)

r = 0;%input("1 for reference; 0 for flight data: "); %reference data: 0 if no, 1 if yes
%lst_data = [Phugoid_V,Phugoid_th,Dutch_yawrate,Dutch_rollrate,D_Dutch_yawrate,D_Dutch_rollrate];
s = input("choose data: ");%3;%input("data source: ");%3;%choose which data set

if r == 0
    lst_dst = [20,10,1,1,1,1]; %minimum distance between peaks
    lst_tres = [1e-4,1e-8,1e-5,1e-5,1e-4,1e-4]; %minimum hight change between peaks
end
if r == 1
    lst_dst = [15,10,1,1,1,1]; %minimum distance between peaks
    lst_tres = [1e-4,1e-8,1e-5,1e-5,3e-4,1e-4]; %minimum hight change between peaks
end

%selects data set
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

%find peaks
t = data(:,1); y = data(:,2);
y = max(max(y)-y,y-min(y));
[PV_y,PV_t] = findpeaks(y,t,'minpeakdistance',dist,'Threshold',Thres);
%remove wrong peaks
if s == 1 && r == 0
    PV_t = PV_t(2:end);
end
if r == 1 && s == 2
    PV_t = PV_t(2:end);
end
if s == 3 || s == 4 || s == 5 || s == 6
    PV_t = PV_t(3:end);
end
if r == 1 && s == 4
    PV_t = PV_t(2:end);
end
%combine peaks and time location
peaks = [];
for q = 1:length(PV_t)
    t_i = PV_t(q);
    num = find(data==t_i);
    x = data(num,2);
    peaks = [peaks ; x];
end
PV_t = [PV_t, peaks];
final =  mean(PV_t(:,2));

%find period
Periods = [];
for q = 1:length(PV_t)-1
    Per  = 2*(PV_t(q+1,1)-PV_t(q,1));
    Periods = [Periods ; Per];
    if q ~= length(PV_t)-1
        Per = (PV_t(length(PV_t),1)-PV_t(q,1))/((-q+length(PV_t))/2);
        Periods = [Periods ; Per];
    end
end
Q = mean(Periods)

%set dimension and velocity
c = 2.0569;
if s == 1 || s == 2
    V0 = TAS_data_P(1,2);
end
if s == 3 || s == 4
    V0 = TAS_data_P(1,2);
end
if s == 5 || s == 6
    V0 = TAS_data_P(1,2);
end

%change to span if assymetric
if s == 3 || s == 4 || s == 5 || s == 6
    final = 0;% data(1,2);
    c = 15.911;
end


%get absolute value around the average
dev_avg = abs(data(:,2)-final);
peaks_avg = abs(PV_t(:,2)-final);

%exponential function dataset from excel. Function: y=ae^(bx)
eq = [];
for q = data(:,1)
    if r == 0 %if flight data
        if s == 1
            a = 7.0437;
            b = -0.004;
            v = 5.5;
            h = 10*Q/16;
            delta = 28;
        end    
        if s == 2
            a = 0.11;
            b = -0.004;
            v = 0.083;
            h = Q/16;
            delta = 16.8;
        end
        if s == 3
            a = 0.2386;
            b = -0.209;
            v = 0.55;
            h = Q/16;
            delta = 5.6;
        end    
        if s == 4
            a = 0.1887;%0.5417
            b = -0.187;%-0.171
            v = 1;
            h = Q/16;
            delta = 6.6;
        end    
        if s == 5
            a = 0.2359;
            b = -0.306;
            v = 1;
            h = Q/16;
            delta = 5;
        end    
        if s == 6
            a = 0.1971;%1.2048
            b = -0.463;%-0.343
            v = 1;
            h = Q/16;
            delta = 6;
        end    
    end
    
    if r == 1 %if reference data
        if s == 1
            a = 7.63;
            b = -0.005;
            v = 5.5;
            h = 10*Q/16;
            delta = 19.3;
        end    
        if s == 2
            a = 0.1083;
            b = -0.005;
            v = 0.083;
            h = 10*Q/16;
            delta = 30;
        end
        if s == 3
            a = 0.3056;
            b = -0.193;
            v = 0.55;
            h = 9.5*Q/16;
            delta = 4.9;
        end    
        if s == 4
            a = 0.2433;%0.5417
            b = -0.208;%-0.171
            v = 1;
            h = 9*Q/16;
            delta = 5.7;
        end    
        if s == 5
            a = 0.262;
            b = -0.264;
            v = 1;
            h = 9*Q/16;
            delta = 4.1;
        end    
        if s == 6
            a = 0.2107;%1.2048
            b = -0.432;%-0.343
            v = 1;
            h = 9*Q/16;
            delta = 5;
        end
    end
    eq = a*exp(b*(q));
    eq_d = a*exp(b*(q-delta));
end

%halftime
x = [data(:,1),eq];
half = log((x(1,2)/2)/a)/b;


%get eigenvalue
V = V0;
eta = 2*pi/Q * c/V;
xi = log(1/2)/half * c/V;

eigenvalue1 = complex(xi,eta);
eigenvalue2 = conj(eigenvalue1);

%using the eigenvalue to replicate the function
test_lst = [];
for p = 1:length(data(:,1))
    p = p/10 + h - delta;
    test = 0.769*a*exp((xi*p)*V/c)*(cos((eta*p)*V/c)+sin((eta*p)*V/c));
    
    test_lst = [test_lst;test];
end


%plotting
figure();
title("peaks")
plot(data(:,1),data(:,2));
hold on
scatter(PV_t(:,1),PV_t(:,2),[],'r');
hold on
plot(data(:,1),test_lst+final);
hold off


figure();
title("amplitudes")
plot(data(:,1),dev_avg);
hold on
scatter(PV_t(:,1),abs(PV_t(:,2)-final),[],'r');
hold on
plot(data(:,1),eq_d);
hold on
plot(data(:,1),abs(test_lst));
hold off