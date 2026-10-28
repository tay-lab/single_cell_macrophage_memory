% Model parameters (Methods, Table 1). Used by TNF_Model_2.
% Source: NFkB Project Andrew\Simulation\20260414 (Apr 2026); copied 2026-10-18 without code changes.
function [knin, Nt, Ki, klin, Kn, kt,gamma,ktl,a,ka20,IKKt,Tdeg,Kca,kb,kd,ka,ki,kp,n1,n2,kdeg,ik_kr,ik_ks,kb3,kb4] = Params

% Tdeg=[2*10^-4,7*10^-4,7.7*10^-4,8.3*10^-4]; % TNF loss  
Tdeg=2*10^-4; % TNF loss     

knin = .01; %1/s Tay
klin = .05; %1/s Tay
Nt = 1; %uM Jensen
Ki = 0.035; %uM Jensen
Kn = 0.029; %uM Jensen
kt = 1.4*10^-7; %1/uM 1/s Ashall
gamma = 0.0063;%1/s Fitted
ktl = 0.5; %1/s Tay
a = 0.0175; %1/uM 1/s Jensen
ka20 = 0.0018; %uM Ashall
IKKt = 2; %uM Jensen
Kca=1;      %uM Assumed

kb=.1;      %Assumed
kd=0.001;   %Assumed

ka = 0.004; %1/min  Ashall
ki = 0.003; %1/min  Ashall
kp = 0.0006; %1/min  Ashall

kdeg=0.00012; %1/s Fitted

ik_kr=0.00002;  %1/s Tay
ik_ks=0.01;     %1/s Tay

n1=2.5; %Assumed
n2=2; %Assumed

kb3=0.001;
kb4=0.0001;

tlIFN=2*10^-5;      %sc IFN translation
tlTnf=1*10^-5;      %sc TNF translation

% knin = 9.4/60; %1/min
% klin = 0.0018/60; %1/min
% Nt = 1; %uM
% Ki = 0.035; %uM
% Kn = 0.029; %uM
% kt = 1.03/60; %1/uM 1/min
% gamma = 0.017/60; %1/min
% ktl = 0.024/60; %1/min
% a = 1.05/60; %1/uM 1/min
% ka = 0.24/60; %1/min
% ki = 0.18/60; %1/min
% kp = 0.036/60; %1/min
% ka20 = 0.0018; %uM
% A20 = 0.0026; %uM
% IKKt = 2.0; %uM


end