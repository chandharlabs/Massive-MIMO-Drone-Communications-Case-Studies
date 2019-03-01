clc;
% close all;
clear all;



%%%%%%%%%%%%% Figure 2a %%%%%%%%%%%%%%%%%%%%%

K = 10:10:1000;
M = 100;
rhop = 100;
B = 20e6;
Bc = 3e6;
fc = 2.4e9;
c = 3e8;
v = 30;
Lambda = (9/10)-(K*2*v*fc/Bc/c);

rhou = 10;
R = B*K.*Lambda.*log2(1+(M*rhou./((rhou*(K-1))+1+((1/(rhou*rhop))*(1+K*rhou)))));
Q = R/1e9;
figure(2);plot(K,Q,'r','Linewidth',2);hold on;

rhou = 1;
R = B*K.*Lambda.*log2(1+(M*rhou./((rhou*(K-1))+1+((1/(rhou*rhop))*(1+K*rhou)))));
Q = R/1e9;
figure(2);plot(K,Q,'g','Linewidth',2);hold on;

rhou = .1;
R = B*K.*Lambda.*log2(1+(M*rhou./((rhou*(K-1))+1+((1/(rhou*rhop))*(1+K*rhou)))));
Q = R/1e9;
figure(2);plot(K,Q,'b','Linewidth',2);

axis([0 1000 1 2.2]);
xlabel('Number of drones');ylabel('Sum Throughput (Gbps)');
legend(3,'Data SNR = 10 dB','Data SNR = 0 dB','Data SNR = -10 dB');


%%%%%%%%%%%%% Figure 3a %%%%%%%%%%%%%%%%%%%%%

M = 100;
rhop = 100;
B = 20e6;
Bc = 3e6;
c = 3e8;
v = 30;

rhoudB = -20:20;
rhou = 10.^(rhoudB/10);

K = 20;
fc = 2.4e9;
Lambda = (9/10)-(K*2*v*fc/Bc/c);
R_K20 = B*K.*Lambda.*log2(1+(M*rhou./((rhou*(K-1))+1+((1./(rhou*rhop)).*(1+K*rhou)))));
Q_K20 = R_K20/1e9;
figure(21);plot(10*log10(rhou),Q_K20);hold on;

fc = 60e9;
Lambda = (9/10)-(K*2*v*fc/Bc/c);
R_K2 = B*K.*Lambda.*log2(1+(M*rhou./((rhou*(K-1))+1+((1./(rhou*rhop)).*(1+K*rhou)))));
Q_K2 = R_K2/1e9;
figure(21);plot(10*log10(rhou),Q_K2,'r--');hold on;
legend('2.4 GHz Carrier Frequency','60 GHz Carrier Frequency');


K = 2;
fc = 2.4e9;
Lambda = (9/10)-(K*2*v*fc/Bc/c);
R_K20 = B*K.*Lambda.*log2(1+(M*rhou./((rhou*(K-1))+1+0*((1./(rhou*rhop)).*(1+K*rhou)))));
Q_K20 = R_K20/1e9;figure(21);plot(rhoudB,Q_K20);hold on;

fc = 60e9;
Lambda = (9/10)-(K*2*v*fc/Bc/c);
R_K2 = B*K.*Lambda.*log2(1+(M*rhou./((rhou*(K-1))+1+0*((1./(rhou*rhop)).*(1+K*rhou)))));
Q_K2 = R_K2/1e9;
figure(21);plot(rhoudB,Q_K2,'r--');
xlabel('Data SNR (dB)');ylabel('Sum Throughput (Gbps)');


%%%%%%%%%%%%% Figure 3b %%%%%%%%%%%%%%%%%%%%%


fc = 60e9;
Pt = .1; % 100 mW
B = [10:500]*1e6;
v = 20;
K = 20;
Bc = 3e6;
M = 145;
N0 = 2e-20; % N0 = 1.38e-23*300*10^(7/10);
rhoudB = -10;
rhou = 10^(rhoudB/10);
Range = sqrt(Pt)*c./(fc*4*pi*sqrt(N0*B*rhou));
Lambda = (9/10)-(K*2*v*fc/Bc/c);
R = K*B*Lambda*log2(1+(M*rhou./((rhou*(K-1))+1+((1./(rhou*rhop)).*(1+K*rhou)))));
Q = R/1e9;
figure(22);[AX,H1,H2] = plotyy(B/1e6,Range,B/1e6,Q);hold on;

rhoudB = 0;
rhou = 10^(rhoudB/10);
Range = sqrt(Pt)*c./(fc*4*pi*sqrt(N0*B*rhou));

figure(22);plot(B/1e6,Range);
set(AX(1),'xlim',[10 500],'ylim',[0 250]);
set(AX(2),'xlim',[10 500],'ylim',[0 40]);hold on;
xlabel('Bandwidth (MHz)'); ylabel(AX(1),'Range (m)');ylabel(AX(2),'Sum Throughput (Mbps)');


%%%%%%%%%%%%% Case Study 1 %%%%%%%%%%%%%%%%%%%%%


% GSD = 20e-2; % m
% FL = 5e-3;
% PS = 2.3e-6/1e3;
% H = GSD*FL/PS;
% OL = .7; 

fc = 2.4e9; % Hz
c = 3e8; % m/s
CR = 200;
NBP = 24;
v = 20;
K = 23;
rhou = 1;
rhop = 100;
Pt = 1; % in W
B = 20e6;
fc = 2.4e9;
rpy = 4096;
rpx = 2048;
FPS = 60;
Q = rpy*rpx*NBP*FPS;
Coherence_Time_Case_Study1 = c/(2*v*fc)/1e-3; % in ms
Bc = 3e6;
Lambda = (9/10)-((2*v*fc*K)/(Bc*c));
Coherence_Interval_Case_Study1_Bc3MHz = Bc*Coherence_Time_Case_Study1*1e-3;
Sum_Throughput_Case_Study1 = K*Q/CR/1e9; % in Gbps
Range = sqrt(Pt)*c./(fc*4*pi*sqrt(N0*B*rhou))/1e3; % in Km
Mreq_Case_Study1_Bc3MHz = (K-1+inv(rhou)+(inv(rhou^2*rhop)*(1+K*rhou)))*((2^(Q/CR/B/Lambda))-1);
 
Bc = 300e3;
Lambda = (9/10)-((2*v*fc*K)/(Bc*c));
Coherence_Interval_Case_Study1_Bc300KHz = Bc*Coherence_Time_Case_Study1*1e-3;
Mreq_Case_Study1_Bc300KHz = (K-1+inv(rhou)+(inv(rhou^2*rhop)*(1+K*rhou)))*((2^(Q/CR/B/Lambda))-1);




%%%%%%%%%%%%% Case Study 2 %%%%%%%%%%%%%%%%%%%%%

fc = 60e9; % Hz
CR = 200;
NBP = 24;
v = 20;
K = 20;
rhou = 1;
rhop = 100;
Pt = 1; % in W
rpy = 4096;
rpx = 2048;
FPS = 60;
Q = rpy*rpx*NBP*FPS*4*4;
Coherence_Time_Case_Study2 = c/(2*v*fc)/1e-3; % in ms
Bc = 3e6;
Lambda = (9/10)-((2*v*fc*K)/(Bc*c));
Coherence_Interval_Case_Study2 = Bc*Coherence_Time_Case_Study2*1e-3;
Sum_Throughput_Case_Study2 = K*Q/CR/1e9; % in Gbps
B = 300e6;
Range_Case_Study2_300MHz = sqrt(Pt)*c./(fc*4*pi*sqrt(N0*B*rhou)); % in Km
Mreq_Case_Study2_300MHz = (K-1+inv(rhou)+(inv(rhou^2*rhop)*(1+K*rhou)))*((2^(Q/CR/B/Lambda))-1);
B = 200e6;
Range_Case_Study2_200MHz = sqrt(Pt)*c./(fc*4*pi*sqrt(N0*B*rhou)); % in Km
Mreq_Case_Study2_200MHz = (K-1+inv(rhou)+(inv(rhou^2*rhop)*(1+K*rhou)))*((2^(Q/CR/B/Lambda))-1);


%%%%%%%%%%%%% Case Study 3 %%%%%%%%%%%%%%%%%%%%%
fc = 5.8e9; % Hz
CR = 1;
NBP = 8;
v = 30;
K = 25;
rhou = 1;
rhop = 100;
Pt = .1; % in W
rpy = 640;
rpx = 480;
FPS = 30;
Q = rpy*rpx*NBP*FPS;
Coherence_Time_Case_Study3 = c/(2*v*fc)/1e-3; % in ms
Bc = 3e6;
Lambda = (9/10)-((2*v*fc*K)/(Bc*c));
Coherence_Interval_Case_Study3 = Bc*Coherence_Time_Case_Study3*1e-3;
B = 20e6;
Sum_Throughput_Case_Study3_K25 = K*Q/CR/1e9; % in Gbps
Range_Case_Study3_K25 = sqrt(Pt)*c./(fc*4*pi*sqrt(N0*B*rhou))/1e3; % in Km
Mreq_Case_Study3_K25 = (K-1+inv(rhou)+(inv(rhou^2*rhop)*(1+K*rhou)))*((2^(Q/CR/B/Lambda))-1);


K = 50;
Sum_Throughput_Case_Study3_K50 = K*Q/CR/1e9; % in Gbps
Range_Case_Study3_K50 = sqrt(Pt)*c./(fc*4*pi*sqrt(N0*B*rhou))/1e3; % in Km
Mreq_Case_Study3_K50 = (K-1+inv(rhou)+(inv(rhou^2*rhop)*(1+K*rhou)))*((2^(Q/CR/B/Lambda))-1);


%%%%%%%%%%%%% Figure 2b %%%%%%%%%%%%%%%%%%%%%

fc = 2.4e9;
lambda = c/fc;
M = 100;
R = 500; % in m
tag_isot = 0;
delta = .5*lambda;

N0 = 2e-20;
B = 20e6;
E_GS = inv(sqrt(2))*[1i ; 1i]; % circularly polarized
E_UAV = inv(sqrt(2))*[1i; 1i]; % circularly polarized
tag_identical = 1; % identically oriented
[z,Pt] = polMismatch_Mag(M,R,tag_identical,tag_isot,delta,E_GS,E_UAV,lambda,N0,B);
figure(22);cdfplot(10*log10(Pt)+30);hold on;
tag_identical = 0;% pseudo-randomly oriented
[z,Pt] = polMismatch_Mag(M,R,tag_identical,tag_isot,delta,E_GS,E_UAV,lambda,N0,B);
figure(22);cdfplot(10*log10(Pt)+30);hold on;


E_GS = [1 ; 0]; % linearly polarized
E_UAV = [1; 0]; % linearly polarized
tag_identical = 1; % identically oriented
[z,Pt] = polMismatch_Mag(M,R,tag_identical,tag_isot,delta,E_GS,E_UAV,lambda,N0,B);
figure(22);cdfplot(10*log10(Pt)+30);hold on;
tag_identical = 0; % pseudo-randomly oriented
[z,Pt] = polMismatch_Mag(M,R,tag_identical,tag_isot,delta,E_GS,E_UAV,lambda,N0,B);
figure(22);cdfplot(10*log10(Pt)+30);hold on;
