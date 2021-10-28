%% TFT example

clc
clear
close all

%% Simulation parameters


    f1 = 60;                % Fundamental frequency
    T = 1/f1;     
    
    Mc = 32;                % Samples per cycle
    fs = f1*Mc;             % Sampling frequency
    Ts = 1/fs;                  

    Nc = 50;                % Number of cycles
    
    n=0:1:Nc*Mc;            % Samples vector
    t=n*Ts;                 % Time vector in sec

%% Test signal

fa = 2;
fc = 2;
c1 = 1; c2 = 0.1;
c3 = 0; c4 = 0.1;

% Values for amp and phase of fundamental component

amp = c1+c2*cos(2*pi*fa*t);
phase= c3+c4*cos(2*pi*fc*t);
% Then frequency 
frec = 60-c4*fc*sin(2*pi*fc*t);

% Signal (you can add harmonic components if desired)
s = (amp.*cos(2*pi*f1*t+phase))';

    
%% TFT estimator for signal s
    harmonic = 1;
    [amp_, phase_, frec_,H0,H1,H2,B,xx] = run_estimate(s,Mc,f1, harmonic);
    
%% Figures

% Frequency responses of the filters related to the 1st harmonic
figure(1)
subplot(3,1,1)
plot(xx,fftshift(abs(fft(H0,4096))),'LineWidth',2)
ylabel('|H_{0,1}(u)|')
grid
subplot(3,1,2)
plot(xx,fftshift(abs(fft(H1,4096))),'LineWidth',2)
ylabel('|H_{1,1}(u)|')
grid
subplot(3,1,3)
plot(xx,fftshift(abs(fft(H2,4096))),'LineWidth',2)
ylabel('|H_{2,1}(u)|')
grid
xlabel('Normalized frequency u = f/f1')

% Estimated values for amplitude, phase and frequency of the 1st harmonic
figure(2)
subplot(4,1,1)
plot(t, s,'LineWidth',2)
axis tight 
grid
ylabel('Input')
subplot(4,1,2)
plot(t,phase,t,phase_,'--','LineWidth',2)
axis([0 t(end) c3-c4 c3+c4])
grid
ylabel('Phase')
legend('Real','Estimated')
subplot(4,1,3)
plot(t,amp,t,amp_,'--','LineWidth',2)
axis tight
axis([0 t(end) c1-c2 c1+c2])
grid
ylabel('Amplitude')
xlabel('Time in sec.')
legend('Real','Estimated')
subplot(4,1,4)
plot(t,frec,t,frec_,'--','LineWidth',2)
axis tight
axis([0 t(end) f1-fc*c4 f1+fc*c4])
grid
ylabel('Frequency')
xlabel('Time in sec.')
legend('Real','Estimated')

%% Functions

function [amp_, phase_, frec_,H0,H1,H2,B,xx] = run_estimate(s,Mc,fc, harmonic)
    %% Parameters
    k=0:1:length(s)-1;
    Ts = fc*Mc;
    scale = (pi)/Ts;
    
    % Create FIR filters
    [est,B,k1,k2,k3,xx] = DTTFs(Mc*4+1,Mc,3,8,2^12);
    
    if mod(Mc,2) == 0 
        offset = 0;
    else
        offset = 1;
    end
    
    % Filters are reversed form of est rows 
    H0 = ((-1)^0)*est(k1((Mc+offset)/2-harmonic),:);     % Order 0
    H1 = ((-1)^1)*est(k2((Mc+offset)/2-harmonic),:);     % Order 1 
    H2 = ((-1)^2)*est(k3((Mc+offset)/2-harmonic),:);     % Order 2
    
    %% Once filter is obtained convolution between the impulse response HX and s signal is done
    Y0z=conv(s,H0);
    Y1z=conv(s,H1);
    Y2z=conv(s,H2);
    Y0=Y0z(1:length(s));
    Y1=Y1z(1:length(s));
    Y2=Y2z(1:length(s));
    
    
   % Then ampitude, phase and frequency is computed
    
    %% Phase needs rotation
    n=0:1:length(Y0)-1;
    ROT=exp(-j*2*pi*harmonic*n'/Mc);
    
    %% Phase and amplitude estimates
    amp_ = 2*abs(Y0);
    phase_ = angle(Y0.*ROT);
    % Frequency estimate
    frec_ =  fc+imag(Y1.*ROT.*exp(-j*phase_))./(scale*amp_);
end

function [est,B,k1,k2,k3,xx] = DTTFs(N,Nf,order,alfa,nmo)
    %% Taylor-Fourier Matrix
    Nh=(N-1)/2; Nfrec=(Nf-1)/2; w1=2*pi/Nf;

    cut=round((order+1)/2);

    B0=zeros(N,cut);
    B0c=zeros(N,order-cut-1);
    BMr=zeros(N,cut*(2*(Nfrec)+1));
    BMi=zeros(N,(order-cut+1)*(2*(Nfrec)+1));
    BM=zeros(N,(2*(Nfrec)+1)*(order+1));
    k=zeros(1,(2*(Nfrec)+1));
    M_Vent=zeros(1,(2*Nh)+1);		
    Vent=zeros(2*Nh+1,2*Nh+1);	
    e=zeros(2*Nh+1,nmo);

    if mod(Nf,2) == 0
    % Taylor-Fourier Matrix for even harmonics
        for kpa = 1:1:(2*(Nfrec)+1)
            for nI=1:1:(cut)
                for n = -Nh : 1: Nh
                    B0((n+Nh+1),nI)=((n^(2*(nI-1)))*exp(j*n*w1*(kpa-1-Nfrec+.5)));
                end
            end
            for nI=1:1:(order+1-cut)
                for n = -Nh : 1: Nh
                    B0c((n+Nh+1),nI)=((n^(2*(nI-1)+1))*exp(j*n*w1*(kpa-1-Nfrec+.5)));
                end    
            end
            BMr(1:N,((kpa-1)*cut+1):(kpa*cut))=B0;
            BMi(1:N,((kpa-1)*(order-cut+1)+1):(kpa*(order-cut+1)))=B0c;

            k(1,kpa)=kpa-1;
        end
        B=[BMr BMi];
    else
        % Taylor-Fourier Matrix for odd harmonics
        for kpa = 1:1:(2*(Nfrec)+1)
            for nI=1:1:(cut)
                for n = -Nh : 1: Nh
                    B0((n+Nh+1),nI)=((n^(2*(nI-1)))*exp(j*n*w1*(kpa-1-Nfrec)));
                end
            end
            for nI=1:1:(order+1-cut)
                for n = -Nh : 1: Nh
                    B0c((n+Nh+1),nI)=((n^(2*(nI-1)+1))*exp(j*n*w1*(kpa-1-Nfrec)));
                end    
            end
            BMr(1:N,((kpa-1)*cut+1):(kpa*cut))=B0;
            BMi(1:N,((kpa-1)*(order-cut+1)+1):(kpa*(order-cut+1)))=B0c;

            k(1,kpa)=kpa-1;
        end
        B=[BMr BMi];
    end
    
    %% Window
    Vent=ventana2d(Nh,alfa);

    %% WLS solution
    est=(inv(B'*Vent*B)*B'*Vent);

    %% Frequency vector
    xx=linspace(-Nf/2,Nf/2, nmo );
   
    %% This indicators are used to determine which rows are related to the different orders
    
    k1=k*cut+1;                                     % order 0
    k2=((2*(Nfrec)+1)*cut)+(k*(order-cut+1))+1;     % order 1
    k3=k*cut+2;                                     % order 2
end

function [V] = ventana2d(Nh,alfa)
    %% The kaiser window is used in the WLS solution
    N=2*Nh+1;
    M_Vent=zeros(1,(2*Nh)+1);		
    V=zeros(2*Nh+1,2*Nh+1);	
    M_Vent=(kaiser(N,alfa))';      
    for nI=1:1:((2*Nh)+1)			
        M_Vent(nI,:)=(kaiser(N,alfa))';     	    
        V(nI,nI)=M_Vent(nI,nI)^(1);   
    end
end
