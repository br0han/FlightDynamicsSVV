run("Parameters.m")
%%%%Data
SDat=load("Elevator.mat"); %Sample Data excluding mass in kg
SMass=load("SampleMasses.mat"); 
SPMass=SMass.PeopleMasses; %Sample Mass of the people in kg


EmptyMass=9165; %(lbs)
EmptyMass=EmptyMass*0.453592; %(kg)
PreFuel=4100; %(lbs)
PreFuel=PreFuel*0.453592; %(kg)

STotPeopleMass=sum(SPMass); %(kg)
STotPreMass=PreFuel+EmptyMass +STotPeopleMass;



SAlt=SDat.SampleData1(:,1);    %Sample Altitude        (ft)
SVelo=SDat.SampleData1(:,2);   %Sample Velocity        (kts)
SAlpha=SDat.SampleData1(:,3);  %Sample Angle of Attack (degrees)
Sde  =SDat.SampleData1(:,4);   % Sample Deflection     (degress)

SFe = SDat.SampleData1(:,7); 
SFFl=SDat.SampleData1(:,7);    %Sample Fuel flow left  (lbs/hr)
SFFr=SDat.SampleData1(:,8);    %Sample Fuel flow right (lbs/hr)
SFUsed=SDat.SampleData1(:,9);  %Sample Fuel used       (lbs)
STemp=SDat.SampleData1(:,10);   %Measured Sample Temperature     (Celsius)

SAlt= 0.3048*SAlt ; %Sample Altitude        (m)
SVelo= 0.514444*SVelo; %Sample Velocity        (m/s)
SAlpha= pi*SAlpha/180; %Sample Angle of Attack (radians)
SFFl= 0.000125998*SFFl; %Sample Fuel flow left  (kg/s)
SFFr= 0.000125998*SFFr; %Sample Fuel flow right (kg/s)
SFUsed= 0.453592*SFUsed; %Sample Fuel used       (kg)
STemp= 273.15+STemp ;  %Measured Sample Temperature     (Kelvin)  



sz=size(SAlt); 

Srho=zeros(sz(1),sz(2)); %Density for Sample data points

SStaticP=zeros(sz(1),sz(2));  %Static pressure Sample data [Pa]

SM=zeros(sz(1),sz(2));      %Mach number Sample data

FM=zeros(sz(1),sz(2));      %Mach number Flight data

SSound=zeros(sz(1),sz(2)); %Speed of sound Sample Data

STempstatic=zeros(sz(1),sz(2)); %Static Air temperature Sample Data

SVeloT=zeros(sz(1),sz(2));     % True speed sample data

SVeloEq=zeros(sz(1),sz(2));     % Equivelent speed sample data



for i=1:sz(1)
    SStaticP(i,1)=p0*(1+(lambda*SAlt(i,1)/Temp0))^(-g/(lambda*R));
    
    SM(i,1) = sqrt(2/(gamma-1)*((1+ p0/SStaticP(i,1) *((1+(gamma-1)*rho0*SVelo(i,1)^2/(2*gamma*p0))^(gamma/(gamma-1)) -1))^((gamma-1)/gamma) -1));
    
    STempstatic(i,1)=STemp(i,1)/(1+(gamma-1)*SM(i,1)^2/2);
    
    SSound(i,1)=sqrt(gamma*R*STempstatic(i,1));
    
    SVeloT(i,1)=SM(i,1)*SSound(i,1);
    
    Srho(i,1)=SStaticP(i,1)/(R* STempstatic(i,1));
    
    SVeloEq(i,1) = SVeloT(i,1) * sqrt(Srho(i,1)/rho0);
    
end

SCL=zeros(sz(1),sz(2)); %Matrix containing CL for Sample data

for i= 1:sz(1)
   SCL(i,1)=2* (STotPreMass-SFUsed(i,1))*g/(Srho(i,1)*SVeloT(i,1)^2*S);
end




delta_cg = 0.2;                                % change in cg  
delta_elev = Sde(7,1) -  Sde(6,1);             % change in elevator deflection 
cm_d = -1/delta_elev * SCL(7:1) *delta_cg/c;




