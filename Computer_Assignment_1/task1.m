clc
clear all

%% Geometry data

a= 11.12; %in m
b=18.7; %in m
c=32.21; %in m
%% Material data
t_s=0.16; %thickness of steel in m
t_a=0.32; %thickness of aluminium in m
t_c=0.40; %thickness of composite in m

ro_s=7850; %density of steel in kg/m^3
ro_a=2700; %density of aluminium kg/m^3
ro_c=1800; %density of composite kg/m^3

E_s=210; %Youngs modulus of steel in GPa
E_a=70; %Youngs modulus of aluminium in GPa
E_c=50; %Youngs modulus of composite in GPa

sigma_y_s=600; %yield stress of steel in MPa
sigma_y_a=500; %yield stress of aluminium in MPa
sigma_y_c=125; %yield stress of composite in MPa

%% Calculation of moment of inetria for steel
y1_s=t_s/2; 
y2_s=(b-(2*t_s))/2+t_s;
y3_s=b-(t_s/2);

A1_s= c*t_s;
A2_s=(b-(2*t_s))*t_s;
A3_s=a*t_s;

c_y_s=((y1_s*A1_s)+(y2_s*A2_s)+(y3_s*A3_s))/(A1_s+A2_s+A3_s);

MI_s=((c*t_s^3)/3)+ (((b-2*t_s)^3)*(2*t_s))/3 + ((t_s^3)*2*a)/3 + ((b-(2*t_s))^2*(t_s*2*a)) + ((c_y_s-t_s)^2*(A1_s+A1_s+A3_s));

%% Calculation of moment of inetria for aluminium
y1_a=t_a/2; 
y2_a=(b-(2*t_a))/2+t_a;
y3_a=b-(t_a/2);

A1_a= c*t_a;
A2_a=(b-(2*t_a))*t_a;
A3_a=a*t_a;

c_y_a=((y1_a*A1_a)+(y2_a*A2_a)+(y3_a*A3_a))/(A1_a+A2_a+A3_a);

MI_a=((c*t_a^3)/3)+ (((b-2*t_a)^3)*(2*t_a))/3 + ((t_a^3)*2*a)/3 + ((b-(2*t_a))^2*(t_a*2*a)) + ((c_y_a-t_a)^2*(A1_a+A1_a+A3_a));

%% Calculation of moment of inetria for composite
y1_c=t_c/2; 
y2_c=(b-(2*t_c))/2 + t_c;
y3_c=b-(t_c/2);

A1_c= c*t_c;
A2_c=(b-(2*t_c))*t_c;
A3_c=a*t_c;

c_y_c=((y1_c*A1_c)+(y2_c*A2_c)+(y3_c*A3_c))/(A1_c+A2_c+A3_c);

MI_c=((c*t_c^3)/3)+ (((b-2*t_c)^3)*(2*t_c))/3 + ((t_c^3)*2*a)/3 + ((b-(2*t_c))^2*(t_c*2*a)) + ((c_y_c-t_c)^2*(A1_c+A1_c+A3_c));