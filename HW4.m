clear all
close all
clc
 
%% Import Input Audio Signal
[Vin,~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice,~]=audioread('outlowsweep.wav');
[OutMidSpice,~]=audioread('outmidsweep.wav');
[OutHighSpice,FsLTSpice]=audioread('outhighsweep.wav');
TsLTSpice=1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact=2;
Fs=FsLTSpice/downSampFact; 

%% Downsample Input Signal
Vin=Vin([1:downSampFact:end]);

%% Sampling Period
Ts=1/Fs;
%% Number of Samples
Nsamp=length(Vin);
%% Simulated time
tstop=Nsamp*Ts;
%% Parameters of Dynamic Element
L1=0.35*10^(-3);
L2=0.35*10^(-3);
L3=3.5*10^(-3);
L4=3.5*10^(-3);
C1=2.8*10^(-6);
C2=2.8*10^(-6);
C3=28*10^(-6);
C4=4.7*10^(-6);
C5=28*10^(-6);
C6=47*10^(-6);
%% Resistive Parameters
R1=10;
RspkLow=8;
R2=10;
RspkMid=8;
RspkHigh=8;

%% WDF setting of free parameters (adaptation conditions) X
% Impedences definition of high frequency
Z6=2*L1/Ts;
Z5=RspkHigh;
Z4=Z6*Z5/(Z6+Z5);
Z3=Z4;
Z2=Ts/(2*C1);
Z1=Z2+Z3;
Z_h=[Z1,Z2,Z3,Z4,Z5,Z6];
% Impedences definition of mid frequency
Z18=R1;
Z17=Ts/(2*C4);
Z16=Z17+Z18;
Z15=Z16;
Z14=RspkMid;
Z13=Z14*Z15/(Z14+Z15);
Z12=Z13;
Z11=2*L3/Ts;
Z10=Z11*Z12/(Z11+Z12);
Z9=Z10;
Z8=Ts/(2*C3);
Z7=Z8+Z9;
Z6=Z7;
Z5=Ts/(2*C2);
Z4=Z5*Z6/(Z5+Z6);
Z3=Z4;
Z2=2*L2/Ts;
Z1=Z2+Z3;
Z_m=[Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10,Z11,Z12,Z13,Z14,Z15,Z16,Z17,Z18];
% Impedences definition of low frequency
Z12=R2; 
Z11=Ts/(2*C6);	
Z10=Z11+Z12;
Z9=Z10;
Z8=RspkLow;
Z7=Z8*Z9/(Z8+Z9);
Z6=Z7;
Z5=Ts/(2*C5);
Z4=Z5*Z6/(Z5+Z6);
Z3=Z4;
Z2=2*L4/Ts;
Z1=Z2+Z3; 
Z_l = [Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10,Z11,Z12]; 

%% Computation of Scattering Matrices X
% Low frequency Scattering Matrices
SS1_l = scattering_series(Z_l(1),Z_l(2),Z_l(3));
SP1_l = scattering_parallel(Z_l(4),Z_l(5),Z_l(6));
SP2_l = scattering_parallel(Z_l(7),Z_l(8),Z_l(9));
SS2_l = scattering_series(Z_l(10),Z_l(11),Z_l(12));
% Mid frequency
SS1_m = scattering_series(Z_m(1),Z_m(2),Z_m(3));
SP1_m = scattering_parallel(Z_m(4),Z_m(5),Z_m(6));
SS2_m = scattering_series(Z_m(7),Z_m(8),Z_m(9));
SP2_m = scattering_parallel(Z_m(10),Z_m(11),Z_m(12));
SP3_m = scattering_parallel(Z_m(13),Z_m(14),Z_m(15));
SS3_m = scattering_series(Z_m(16),Z_m(17),Z_m(18));
% High frequency
SS1_h = scattering_series(Z_h(1),Z_h(2),Z_h(3));
SP1_h = scattering_parallel(Z_h(4),Z_h(5),Z_h(6));

%% Initialization of Waves
% High 
b6_h = 0; %L
b2_h = 0; %C
a5_h = 0; %R
b1_h = 0; %Vin
a1_h = 0;
a4_h = 0;
% Mid
b2_m = 0;
b5_m = 0;
b8_m = 0;
b11_m = 0;
b17_m = 0;
a14_m = 0;
a18_m = 0;

a13_m = 0;
a10_m = 0;
a7_m = 0;
a4_m = 0;
a1_m = 0;
a16_m = 0;
% Low
a12_l = 0; % R
b11_l = 0; % C
a8_l = 0; % R
b5_l = 0; % C
b2_l = 0; % L
b1_l = 0;
a1_l = 0;
a4_l = 0;
a10_l = 0;
a7_l = 0;

%% Initialize Output Signals
% Low
VoutLow=zeros(size(Vin));
% Mid
VoutMid=zeros(size(Vin));
% High
VoutHigh=zeros(size(Vin));

ii=0;
while (ii<Nsamp)
    ii=ii+1;

    %% Manage Dynamic Elements
    % High
    a6_h = -b6_h;
    a2_h = b2_h;
    % Mid
    a2_m = -b2_m;
    a5_m = b5_m;
    a8_m = b8_m;
    a11_m = -b11_m;
    a17_m = b17_m; 
    % Low
    a2_l = -b2_l;
    a5_l = b5_l;
    a11_l = b11_l;

    %% Forward Scan
    % High
    b4_h = SP1_h(1,:) * [a4_h;a5_h;a6_h];
    a3_h = b4_h;
    b1_h = SS1_h(1,:) * [a1_h;a2_h;a3_h];
    % Mid
    b16_m = SS3_m(1,:) * [a16_m;a17_m;a18_m]; 
    a15_m = b16_m;
    b13_m = SP3_m(1,:) * [a13_m;a14_m;a15_m]; 
    a12_m = b13_m;
    b10_m = SP2_m(1,:) * [a10_m;a11_m;a12_m]; 
    a9_m = b10_m;
    b7_m = SS2_m(1,:) * [a7_m;a8_m;a9_m];
    a6_m = b7_m;
    b4_m = SP1_m(1,:) * [a4_m;a5_m;a6_m]; 
    a3_m = b4_m;
    b1_m = SS1_m(1,:) * [a1_m;a2_m;a3_m]; 
    % Low
    b10_l = SS2_l(1,:) * [a10_l;a11_l;a12_l];
    a9_l = b10_l;
    b7_l = SP2_l(1,:) * [a7_l;a8_l;a9_l];
    a6_l = b7_l;
    b4_l = SP1_l(1,:) * [a4_l;a5_l;a6_l];
    a3_l = b4_l;
    b1_l = SS1_l(1,:) * [a1_l;a2_l;a3_l];

    %% Local Root Scattering
    % High
    a1_h = 2*Vin(ii) - b1_h;
    % Mid
    a1_m = 2*Vin(ii) - b1_m;
    % Low
    a1_l = 2*Vin(ii) - b1_l;

    %% Backward Scan
    % High
    b3_h = SS1_h(3,:) * [a1_h;a2_h;a3_h];
    a4_h = b3_h;  
    b5_h = SP1_h(2,:) * [a4_h;a5_h;a6_h];

    b2_h = SS1_h(2,:) * [a1_h;a2_h;a3_h];
    b6_h = SP1_h(3,:) * [a4_h;a5_h;a6_h];

    % Mid
    b3_m = SS1_m(3,:) * [a1_m;a2_m;a3_m];
    b2_m = SS1_m(2,:) * [a1_m;a2_m;a3_m];
    a4_m = b3_m;  
    b6_m = SP1_m(3,:) * [a4_m;a5_m;a6_m];
    b5_m = SP1_m(2,:) * [a4_m;a5_m;a6_m];
    a7_m = b6_m;
    b9_m = SS2_m(3,:) * [a7_m;a8_m;a9_m];
    b8_m = SS2_m(2,:) * [a7_m;a8_m;a9_m];
    a10_m = b9_m;
    b12_m = SP2_m(3,:) * [a10_m;a11_m;a12_m];
    b11_m = SP2_m(2,:) * [a10_m;a11_m;a12_m];
    a13_m = b12_m;
    b15_m = SP3_m(3,:) * [a13_m;a14_m;a15_m];
    b14_m = SP3_m(2,:) * [a13_m;a14_m;a15_m];
    a16_m = b15_m;
    b18_m = SS3_m(3,:) * [a16_m;a17_m;a18_m];
    b17_m = SS3_m(2,:) * [a16_m;a17_m;a18_m];

    % Low
    b3_l = SS1_l(3,:) * [a1_l;a2_l;a3_l];
    a4_l = b3_l;
    b6_l = SP1_l(3,:) * [a4_l;a5_l;a6_l];
    a7_l = b6_l;
    b9_l = SP2_l(3,:) * [a7_l;a8_l;a9_l];
    a10_l = b9_l;
    b12_l = SS2_l(3,:) * [a10_l;a11_l;a12_l];

    b8_l = SP2_l(2,:) * [a7_l;a8_l;a9_l];
    b2_l = SS1_l(2,:) * [a1_l;a2_l;a3_l];
    b5_l = SP1_l(2,:) * [a4_l;a5_l;a6_l];
    b11_l = SS2_l(2,:) * [a10_l;a11_l;a12_l];

    %% Read Output X
    % High
    VoutHigh(ii)=-(a5_h+b5_h)/2;
    % Mid
    VoutMid(ii)=(a14_m+b14_m)/2;
    % Low
    VoutLow(ii)=-(a8_l+b8_l)/2;    
    
end

%% Output Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(TsLTSpice*[1:length(OutLowSpice)],OutLowSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutLow,'b--','Linewidth',1); grid on; xlim([0,tstop]); 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
plot(TsLTSpice*[1:length(OutMidSpice)],OutMidSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutMid,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(TsLTSpice*[1:length(OutHighSpice)],OutHighSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutHigh,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% Error Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(Ts*[1:Nsamp],OutLowSpice([1:downSampFact:end])-VoutLow,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
title(['Error Signals. $F_{\mathrm{s}}=$ ',num2str(Fs),' Hz'],'Fontsize',18,'interpreter','latex');
subplot(312)
plot(Ts*[1:Nsamp],OutMidSpice([1:downSampFact:end])-VoutMid,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(Ts*[1:Nsamp],OutHighSpice([1:downSampFact:end])-VoutHigh,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% compute scattering matrix  X
function S = scattering_series(imp1, imp2, imp3)
   S = eye(3) - (2/(imp1+imp2+imp3))*[imp1; imp2; imp3]*[1,1,1];
end
function S = scattering_parallel(imp1, imp2, imp3)
   S = ((2/(1/imp1+1/imp2+1/imp3))*[1/imp1;1/imp2;1/imp3]*[1,1,1] - eye(3)).';
end 