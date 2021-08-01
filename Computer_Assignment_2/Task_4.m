%% Task_4

Result=zeros(17,6);
 
 
 for i=1:30/2+1

Q=150*10^3*9.81;     %[N/m]
L=28+i*2;                %[m]
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
%Sigma_T=max(Bi)*omega/Kw;

%  Result(i,:)=[M alfa Sigma_T Sigma_M sigma_max Sigma_M/Sigma_T ];
Result(i,:)=[M alfa Sigma_T Sigma_M sigma_max Sigma_M/Sigma_T ];

 end
x = 30:2:60;
y = Result(1:16,1);
plot(x,transpose(y));
xlabel('Length [m]','interpreter','latex')
ylabel('Bending moment [N-m]','interpreter','latex')
title('Change in the Maximum Bending moment','interpreter','latex')

x = 30:2:60;
y = Result(1:16,2);
plot(x,transpose(y));
xlabel('Length [m]','interpreter','latex')
ylabel('\alpha')
title('Change in the Mixed torsion parameter','interpreter','latex')

x = 30:2:60;
y = Result(1:16,4);
plot(x,transpose(y));
xlabel('Length [m]','interpreter','latex')
ylabel('\sigma [MPa]')
title('Change in the Maximum Bending stress','interpreter','latex')

x = 30:2:60;
y = Result(1:16,3);
plot(x,transpose(y));
xlabel('Length [m]','interpreter','latex')
ylabel('\sigma_{w} [MPa]')
title('Change in the Maximum Warping stress','interpreter','latex')

x = 30:2:60;
y = Result(1:16,5);
plot(x,transpose(y));
xlabel('Length [m]','interpreter','latex')
ylabel('\sigma^{max}_{T} [MPa]')
title('Change in the Maximum absolute total normal stress','interpreter','latex')
