
SThP = load('STP.mat');
FThP = load('FTP.mat');

stps = SThP.SampleThrustParams;
ftps = FThP.FlightThrustParams;

%format [hp M FFl FFr DeltaT]

S_Thrust = zeros(size(stps, 1), 2);
F_Thrust = zeros(size(ftps, 1), 2);

for i = 1:size(stps,1)
    tempstp = stps(i, :);
    dlmwrite('matlab.dat', tempstp, 'delimiter', ' ');
    system('FindThrust.exe');
    S_Thrust(i, :) = load('thrust.dat');
end

for j = 1:size(ftps,1)
    tempftp = ftps(i, :)
    dlmwrite('matlab.dat', tempftp, 'delimiter', ' ');
    system('FindThrust.exe');
    F_Thrust(i, :) = load('thrust.dat');
end



