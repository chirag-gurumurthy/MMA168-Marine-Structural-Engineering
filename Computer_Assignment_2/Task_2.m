
clc, clear all, clf

Q=150*10^3*9.81;     %[N/m]
L=45;                %[m]
B=30;                %[m]
H=12;                %[m]
C=6;                 %[m]
D=4;                 %[m]
t=25*10^-3;          %[m]
e=4.36;              %[m]
Kv=1.031*10^-3;      %[m^4]
Iy=40.58;            %[m]
%material properties
E=210*10^9;
v=.3;
G=E/(2*(v+1));
%cross section properties
beta=H/B;
gama=C/B;
epsilon=beta*(2*gama*(3-4*gama^2)+3*beta)/(1+2*gama*(4*gama^2-6*gama+3)+6*beta);
omega1=epsilon*B^2/2;
omega2=(beta-epsilon)*B^2/2;
omega3=omega2+gama*(beta+epsilon)*B^2;
Kw=(B*t/3)*((1+2*epsilon)*omega1^2+2*(beta+gama-epsilon)*omega2^2+2*gama*omega3*(omega2+omega3));
%vlasov or st: venant
alfa=abs(-G*Kv*L^2/(pi^2*E*Kw));
%bi moment
B1=Q*(B/2-D)*L^2/12;
B3=Q*(B/2-D)*L^2/24;
Bi=[B1 B3];
omega=[omega1 omega2 omega3];
%normal stress from torsion
Sigma_T=max(Bi)*max(omega)/Kw;
%bending moment from load
M=-Q*L^2/12;
%stress from bending
if abs(-e)>=abs(H-e)
    Sigma_M=M/Iy*-e;
elseif abs(H-e)>=abs(-e)
    Sigma_M=M/Iy*(H-e);
end

sigma_max=-Sigma_T+Sigma_M;
Sigma_T=max(Bi)*omega/Kw;

z1=linspace(-e ,0);
z2=linspace(0 ,H-e);
sigma_m1= M/Iy*z1*10^-6;
sigma_m2= M/Iy*z2*10^-6;
scale4=10^-7;

scale=1.2;
figure(1)
Y=[C, 0, 0, B, B, B-C];
Z=[H, H, 0, 0 H, H,];
plot(Y,Z,'black-','LineWidth',3)
axis([-B*(scale-1)-1 B*scale+1 -H*(scale-1)-1 H*scale+1])
hold on
plot([0],[e],'oG','LineWidth',6)
plot([B],[e],'oG','LineWidth',6)
plot(sigma_m1*.3+B/2,z1+e,'r')
plot(sigma_m2*.3+B/2,z2+e,'b')
plot([0,B],[e e],'--black')
A=[num2str(min(sigma_m2)),' [MPa]']   
text(0,H+1,A);
B1=[num2str(max(sigma_m1)),' [MPa]']   
text(B-1,0-1,B1);
xlabel('Breadth [m]','interpreter','latex')
ylabel('Height [m]','interpreter','latex')
title('Bending moment load','interpreter','latex')
legend('Hull cross-section','Zero stress point','Neutral axis','Tensile stress','Compressive stress')
hold on

figure(2)
Y=[C, 0, 0, B, B, B-C];
Z=[H, H, 0, 0 H, H,];
plot(Y,Z,'black-','LineWidth',3)
axis([-B*(scale-1)-1 B*scale+1 -H*(scale-1)-1 H*scale+6   ])
hold on
scale2=1*10^-7;

plot([0,B/2],[-Sigma_T(1)*scale2,0],'r')
plot([B/2,B],[0,Sigma_T(1)*scale2],'b')
plot([0,(Sigma_T(1)+Sigma_T(2))*scale2],[e,H],'b')
plot([-Sigma_T(1)*scale2,0],[0,e],'r')
plot([B-(Sigma_T(1))*scale2,B+0],[0,e],'b')
plot([B-0,B+(Sigma_T(1)+Sigma_T(2))*scale2],[e,H],'r')
plot([0, C],[H+(Sigma_T(3)*0-Sigma_T(2))*scale2, H-(Sigma_T(3))*scale2],'b')
plot([B, B-C],[H+(Sigma_T(3)*0+Sigma_T(2))*scale2, H+(Sigma_T(3))*scale2],'r')

 T1=[num2str((-Sigma_T(3))*10^-6),' [MPa]'];   
 text(C,H-(Sigma_T(3))*scale2,T1);
 T2=[num2str(Sigma_T(3)*10^-6),' [MPa]'];
 text(B-C+1,H+(Sigma_T(3))*scale2,T2);


plot([B/2],[0],'oG','LineWidth',6)
plot([0],[e],'oG','LineWidth',6)
plot([B],[e],'oG','LineWidth',6)
xlabel('Breadth [m]','interpreter','latex')
ylabel('Height [m]','interpreter','latex')
title('Torsion load distrubution','interpreter','latex')
legend('Hull cross-section','Zero stress point','Neutral axis','Tensile stress','Compressive stress')


figure(3)
sigma_m1= M/Iy*z1;
sigma_m2= M/Iy*z2;
Y=[C, 0, 0, B, B, B-C];
Z=[H, H, 0, 0 H, H,];
plot(Y,Z,'black-','LineWidth',3)
%axis([-B*(scale-1)-1 B*scale+1 -H*(scale-1)-2 H*scale+3   ])
grid on
grid minor
hold on
scale2=1*10^-7;


plot([0,B],[(-Sigma_T(1)+min(sigma_m2))*scale2,(Sigma_T(1)+min(sigma_m2))*scale2],'-m')
plot([(-Sigma_T(1)+min(sigma_m2))*scale2,(Sigma_T(1)+Sigma_T(2)+max(sigma_m1))*scale2],[0,H],'-m')
plot([B+(-Sigma_T(1)+min(sigma_m2))*scale2,   B+(Sigma_T(1)+Sigma_T(2)+max(sigma_m1))*scale2],[0,H],'-m')
plot([0, C],[H+(-Sigma_T(2)+min(sigma_m2))*scale2, H+(-Sigma_T(3)+min(sigma_m2))*scale2],'-m')
plot([B, B-C],[H+(Sigma_T(2)+min(sigma_m2))*scale2, H+(Sigma_T(3)+min(sigma_m2))*scale2],'-m')

A1=[num2str(((-Sigma_T(3)+min(sigma_m2)))*10^-6),' [MPa]'];   
 text(C+1,H+(-Sigma_T(3)+min(sigma_m2))*scale2,A1);
 plot(C,H+(-Sigma_T(3)+min(sigma_m2))*scale2,'*r','LineWidth',3)
 plot(C,H,'*r','LineWidth',3)
 
xlabel('Breadth [m]','interpreter','latex')
ylabel('Height [m]','i1nterpreter','latex')
title('Combined load case','interpreter','latex')






