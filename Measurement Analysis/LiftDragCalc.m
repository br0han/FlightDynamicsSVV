run("Parameters.m")
%%%%Data
SDat=load("SampleData.mat"); %Sample Data excluding mass in kg
SMass=load("SampleMasses.mat"); 
SPMass=SMass.PeopleMasses; %Sample Mass of the people in kg


FDat=load("LiftFlightData.mat"); %Measurements from flight
FMass=load("FlightMassData"); %Mass data from flight
FPMass=FMass.FlightData;

EmptyMass=9165; %(lbs)
EmptyMass=EmptyMass*0.453592; %(kg)
PreFuel=4100; %(lbs)
PreFuel=PreFuel*0.453592; %(kg)


STotPeopleMass=sum(SPMass); %(kg)
FTotPeopleMass=sum(FPMass);

STotPreMass=PreFuel+EmptyMass +STotPeopleMass;
FTotPreMass=PreFuel+EmptyMass +FTotPeopleMass;

%SAMPLE DATA
SAlt=SDat.SampleData(:,1);    %Sample Altitude        (ft)
SVelo=SDat.SampleData(:,2);   %Sample Velocity        (kts)
SAlpha=SDat.SampleData(:,3);  %Sample Angle of Attack (degrees)
SFFl=SDat.SampleData(:,4);    %Sample Fuel flow left  (lbs/hr)
SFFr=SDat.SampleData(:,5);    %Sample Fuel flow right (lbs/hr)
SFUsed=SDat.SampleData(:,6);  %Sample Fuel used       (lbs)
STemp=SDat.SampleData(:,7);   %Measured Sample Temperature     (Celsius)

SAlt= 0.3048*SAlt ; %Sample Altitude        (m)
SVelo= 0.514444*SVelo; %Sample Velocity        (m/s)
SAlpha= pi*SAlpha/180; %Sample Angle of Attack (radians)
SFFl= 0.000125998*SFFl; %Sample Fuel flow left  (kg/s)
SFFr= 0.000125998*SFFr; %Sample Fuel flow right (kg/s)
SFUsed=0.453592*SFUsed; %Sample Fuel used       (kg)
STemp=273.15+STemp ;  %Measured Sample Temperature     (Kelvin)  

%FLIGHT DATA
FAlt=FDat.LiftFlightData(:,1);  %Flight Altitude (m)
FVelo=FDat.LiftFlightData(:,2); %Flight Velocity (m/s)
FAlpha=FDat.LiftFlightData(:,3); %Flight Angle of Attack (degrees)
FFFl=FDat.LiftFlightData(:,4);   %Flight Fuel Flow left (lbs/hr)
FFFr=FDat.LiftFlightData(:,5);   %Flight Fuel flow right (lbs/hr)
FFUsed=FDat.LiftFlightData(:,6);  %Flight Fuel Used (lbs)
FTemp=FDat.LiftFlightData(:,7);   %Measured Flight Temperature (Celsius)

FAlt=0.3048*FAlt ; %Flight Altitude (m)
FVelo=0.514444*FVelo ; %Flight Velocity (m/s)
FAlpha=pi*FAlpha/180 ; %Flight AoA (radians)
FFFl=0.000125998*FFFl ; %Flight Fuel Flow Left (kg/s)
FFFr=0.000125998*FFFr ; %Flight Fuel Flow right (kg/s)
FFUsed=0.453592*FFUsed ; %Flight Fuel  used (kg)
FTemp=FTemp+273.15 ; %Measured Flight Temperature (Kelvin)


STempISA=[278.23; 278.21; 278.21; 278.21; 278.21; 278.04];
FTempISA=[264.4; 264.42; 264.4; 264.4; 264.4; 264.4;];

sz=size(SAlt); 

Srho=zeros(sz(1),sz(2)); %Density for Sample data points
Frho=zeros(sz(1),sz(2)); %Density for Flight Data points

SDeltaTemp= STemp-STempISA;
FDeltaTemp= FTemp-FTempISA;

SStaticP=zeros(sz(1),sz(2));  %Static pressure Sample data [Pa]
FStaticP=zeros(sz(1),sz(2));  %Static Pressure flight data [Pa]

SM=zeros(sz(1),sz(2));      %Mach number Sample data
FM=zeros(sz(1),sz(2));      %Mach number Flight data

SSound=zeros(sz(1),sz(2)); %Speed of sound Sample Data
FSound=zeros(sz(1),sz(2)); %Speed of sound flight data

STempstatic=zeros(sz(1),sz(2)); %Static Air temperature Sample Data
FTempstatic=zeros(sz(1),sz(2)); %Static Air temperature Flight Data

SVeloT=zeros(sz(1),sz(2));
FVeloT=zeros(sz(1),sz(2));

for i=1:sz(1)
    SStaticP(i,1)=p0*(1+(lambda*SAlt(i,1)/Temp0))^(-g/(lambda*R));
    FStaticP(i,1)=p0*(1+(lambda*FAlt(i,1)/Temp0))^(-g/(lambda*R));

    SM(i,1) = sqrt(2/(gamma-1)*((1+ p0/SStaticP(i,1) *((1+(gamma-1)*rho0*SVelo(i,1)^2/(2*gamma*p0))^(gamma/(gamma-1)) -1))^((gamma-1)/gamma) -1));
    FM(i,1) = sqrt(2/(gamma-1)*((1+ p0/FStaticP(i,1) *((1+(gamma-1)*rho0*FVelo(i,1)^2/(2*gamma*p0))^(gamma/(gamma-1)) -1))^((gamma-1)/gamma) -1));
    
    STempstatic(i,1)=STemp(i,1)/(1+(gamma-1)*SM(i,1)^2/2);
    FTempstatic(i,1)=FTemp(i,1)/(1+(gamma-1)*FM(i,1)^2/2);
    
    SSound(i,1)=sqrt(gamma*R*STempstatic(i,1));
    FSound(i,1)=sqrt(gamma*R*FTempstatic(i,1));
    
    SVeloT(i,1)=SM(i,1)*SSound(i,1);
    FVeloT(i,1)=FM(i,1)*FSound(i,1);
    
    Srho(i,1)=SStaticP(i,1)/(R* STempstatic(i,1));
    Frho(i,1)=FStaticP(i,1)/(R* FTempstatic(i,1));
end


FThrustL =[3346.92; 2480.35; 2233.03 ;1917.33; 2125.58; 1924.03 ]; %Left Engine Thrust for Flight data
FThrustR =[3684.7; 2836.5; 2517.5; 2116.72; 2432.98; 2223.14];   %Right Engine Thrust for Flight data

FThrust = FThrustL+FThrustR;

SCL=zeros(sz(1),sz(2)); %Matrix containing CL for Sample data
FCL=zeros(sz(1),sz(2)); %Matrix containing CL for Flight data
for i= 1:sz(1)
   SCL(i,1)=2* (STotPreMass-SFUsed(i,1))*g/(Srho(i,1)*SVeloT(i,1)^2*S);
   FCL(i,1)=2* (FTotPreMass-FFUsed(i,1))*g/(Frho(i,1)*FVeloT(i,1)^2*S);
end

FCD=zeros(sz(1),sz(2));
Glider=zeros(sz(1),sz(2));
for i=1:sz(1)
   FCD(i,1)= 2*FThrust(i,1)/(Frho(i,1)*FVeloT(i,1)^2*S);
   Glider(i,1)=FCL(i,1)/FCD(i,1);
end


FCL2=FCL.^2

FCLMatrix=ones(sz(1),2);
FCLMatrix(1:sz(1),1)=FAlpha;

FCLSol=inv(FCLMatrix' *FCLMatrix)*FCLMatrix'*FCL;
FClalpha=FCLSol(1,1);
Clalpha0=FCLSol(2,1);


FCDMatrix=ones(sz(1),2);
FCDMatrix(1:sz(1),1)=FCL2;

FCDSol=inv(FCDMatrix' *FCDMatrix)*FCDMatrix'*FCD;
CD0=FCDSol(2,1);
e=1/(pi*FCDSol(1,1)*b^2/S)

plot(FCD,FAlpha)
