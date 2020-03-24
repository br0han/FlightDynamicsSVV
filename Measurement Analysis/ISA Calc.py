from math import *

print("International Standard Atmosphere Calculator \n") 

#Values of Standard constants
g=9.80665
R=287



#Choosing the input unit
print("1. Calculate ISA for meters")
print("2. Calculate ISA for feet")
print("3. Calculate ISA for FL")

choice=int(input("Enter your choice "))

#Converting the input to meters
if choice is 1:
    alt=float(input("At which altitude would you like to know the properties of air? (in meters) "))

if choice is 2:
    height=float(input("At which altitude would you like to know the properties of air? (in feet) "))
    alt=0.3048*height

if choice is 3:
    height=float(input("At which altitude would you like to know the properties of air? (in FL) "))
    alt=30.48*height

#Choosing the sea level conditions 
print("1. Use standard sea-level conditions ")
print("2. Use different sea-level conditions ")
cond=int(input("Enter your choice "))

if cond is 1:
    p0= 101325
    T0= 288.16

if cond is 2:
    p=float(input("Enter your sea-level Pressure "))
    T=float(input("Enter your sea-level Temperature "))

p=p0
T=T0
y=0
dy=0.1


def Gradient(P,t,A,DY):
    P=P*(((t+(A*DY))/t)**(-g/(A*R)))
    t=t+(A*DY)
    return [t,P]
    #Gradient region

def Isothermal(P,t,DY):
    t=t
    P=P*exp((-g*DY)/(R*t))
    return [t,P]
    #Isothermal region

while not y>=alt:
    y=y+dy
    if 0<y<11000:
        a=-0.0065
        T=Gradient(p,T,a,dy)[0]
        p=Gradient(p,T,a,dy)[1]
    elif 11000<y<20000:
        T=Isothermal(p,T,dy)[0]
        p=Isothermal(p,T,dy)[1]
    elif 20000<y<32000:
        a=0.001
        T=Gradient(p,T,a,dy)[0]
        p=Gradient(p,T,a,dy)[1]
    elif 32000<y<47000:
        a=0.0028
        T=Gradient(p,T,a,dy)[0]
        p=Gradient(p,T,a,dy)[1]
    elif 47000<y<51000:
        T=Isothermal(p,T,dy)[0]
        p=Isothermal(p,T,dy)[1]
    elif 51000<y<71000:
        a=-0.0028
        T=Gradient(p,T,a,dy)[0]
        p=Gradient(p,T,a,dy)[1]
    elif 71000<y<86000:
        a=-0.002
        T=Gradient(p,T,a,dy)[0]
        p=Gradient(p,T,a,dy)[1]

if alt>86000:
    print("Sorry, but the International Standard Atmosphere only computes values up to 86000m")
elif alt<=86000:
    rho=p/(R*T)
    print()
    print("Altitude: ", alt ,"m")
    print("Temperature: ", round(T,2), "K")
    print("Pressure: ", round(p,3), "Pa", "(",round(100*p/p0),"% SL)")
    print("Density: ", round(rho,6), "kg/(m^3)")
    input()

print("\nReady")
