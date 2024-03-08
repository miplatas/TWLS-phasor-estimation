%% msk_v1
% MSK como CPFSK con h=1/2

clc
clear
close all

%% Parametros

    % Mensaje a transmitir
    I=[1 -1 1 -1 -1 -1 1 1 1 1 1 -1 1 -1 1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 -1 1 1 -1 -1 -1]'; 



    fc = 2400000000;               % Frecuencia de la portadora (2.4 GHz)
    Tc = 1/fc;

    f = fc/218;                 % Frecuencia de transmision de simbolos (11 Mbp - 218, 54 Mbps - 45)
    T = 1/f;     
    
    Mc = 100;
    fs = fc*Mc;            % Frecuencia de muestreo (100 muestras por ciclo)
    Ts = 1/fs;                  


%% Modulacion MSK

    % Tipo de pulso (1 RC, 0 Rect)
    p = 0;  
    
    % Modulacion
    [s, d, fase, t] = msk(I,p,f,T,fc,Tc,fs,Ts);

    
%% Representacion polinomial (un tono, luego multitono)
    harmonic = 1;
    [amp_, phase_, frec_,H0,H1,H2,B,xx] = run_estimate(s,Mc,fc, harmonic);
    
%% Figuras

% Figura 1
close all
figure(1)
title('Modulador')
subplot(3,1,1)
plot(t, d,'LineWidth',2)
axis tight 
grid
ylabel('Mensaje')
subplot(3,1,2)
plot(t,fase,'LineWidth',2)
axis tight
grid
ylabel('Fase')
subplot(3,1,3)
plot(t, s,'LineWidth',2)
axis tight 
grid
ylabel('Modulada')
xlabel('Tiempo en seg.')

figure(2)
subplot(3,1,1)
plot(t, s,'LineWidth',2)
axis tight 
grid
ylabel('Modulada')
subplot(3,1,2)
plot(t,fase,t,phase_,'--','LineWidth',2)
axis tight
grid
ylabel('Fase')
legend('Teorica','Estimada')
subplot(3,1,3)
plot(t,fc+d*(1/(4*T)),t,frec_,'--','LineWidth',2)
axis([0 max(t) fc-2.1*1/(4*T) fc+2.1*1/(4*T)])
grid
ylabel('Frecuencia')
xlabel('Tiempo en seg.')
legend('Teorica','Estimada')

%% Funciones usadas

function [s, d, f, t1D] = msk(I,p,f,T,fc,Tc,fs,Ts)
%% Modulacion MSK de mensaje

% I -> vector con elementos 1, -1 que contiene informacion a transmitir
% p -> Selector de pulso a usar
% s <- Señal modulada
% d <- Pulsos a transmitir
% f <- Fase usada para modular s
% t1D <- vector de tiempo
% Mc <- Muestras por ciclo de portadora
% Ts <- Tiempo de muestreo

    %% Parametros de modulacion

    Ns = numel(I);          % Numero de simbolos
    Mc = Tc/Ts;             % Muestras por ciclo de portadora
    Cs = T/Tc;              % Ciclos de portadora por simbolo
    Ms = Mc*Cs;             % Muestras por simbolo
    k = 0:1:Ms-1;           % Tiempo pulso muestras
    t = k*Ts;               % Tiempo pulso seg
    N = Ns*Ms;              % Muestras totales simulacion
    n = 0:1:N-1;            % Tiempo simulacion muestras
    t1D = n*Ts;             % Tiempo simulacion seg
    
    %% Seleccion de tipo de pulso a usar
    
    switch p 
        case 1
            % Pulso Raise Cosine
            g = (1)/(2*T)*(1-cos(2*pi*t/T));    
            q = (t-T*sin(2*pi*t/T)/(2*pi))/(2*T);
        otherwise
            % Pulso rectangular
            g = cos(2*pi*f*t*0);    
            q = 0.5*t/max(t);
    end

    %% Modulacion
    
    t2D = (reshape(t1D,[round(Ms), round(Ns)]))';

    Memoria = 0;
    for r = 1:1:Ns
        d2(r,:) = I(r)*g; 
        fase2(r,:) = 0.5*pi*Memoria + pi*I(r)*q;
        Memoria = Memoria + I(r);
        s2(r,:) = cos(2*pi*fc*t2D(r,:) + fase2(r,:));
    end

    %% Salidas
    
    % Señal modulada
    s=reshape(s2',[numel(s2), 1]); 
    
    % Tren de pulsos
    d=reshape(d2',[numel(d2), 1]);
    
    % Fase usada para modular
    f=reshape(fase2',[numel(fase2), 1]);
end

function [amp_, phase_, frec_,H0,H1,H2,B,xx] = run_estimate(s,Mc,fc, harmonic)
    k=0:1:length(s)-1;
    [est,B,k1,k2,k3,xx] = DTTFs(Mc*4+1,Mc,3,8,2^12);
    H0 = est(k1(Mc/2-harmonic),:);     % Filtro para estimación 1 arm, orden 0
    H1 = est(k2(Mc/2-harmonic),:);     % Filtro para estimación 1 arm, orden 1
    H2 = est(k3(Mc/2-harmonic),:);     % Filtro para estimación 1 arm, orden 2
    
    Y0z=conv(H0,s);
    Y1z=conv(s,H1);
    Y2z=conv(s,H2);
    Y0=Y0z(1:length(s));
    Y1=Y1z(1:length(s));
    Y2=Y2z(1:length(s));
    
    n=0:1:length(Y0)-1;
    ROT=exp(-j*2*pi*harmonic*n'/Mc);
    
    amp_ = 2*abs(Y0);
    phase_ = angle(Y0.*ROT);
    scale = -1.3091e-11;
    frec_ =  fc+imag(Y1.*ROT.*exp(-j*phase_))./(scale*amp_);
end

function [est,B,k1,k2,k3,xx] = DTTFs(N,Nf,order,alfa,nmo)

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
    % Se genera la matriz de Taylor-Fourier
        for kpa = 1:1:(2*(Nfrec)+1)
            for nI=1:1:(cut)
                for n = -Nh : 1: Nh
                    B0((n+Nh+1),nI)=double((n^(2*(nI-1)))*exp(j*n*w1*(kpa-1-Nfrec+.5)));
                end
            end
            for nI=1:1:(order+1-cut)
                for n = -Nh : 1: Nh
                    B0c((n+Nh+1),nI)=double((n^(2*(nI-1)+1))*exp(j*n*w1*(kpa-1-Nfrec+.5)));
                end    
            end
            BMr(1:N,((kpa-1)*cut+1):(kpa*cut))=B0;
            BMi(1:N,((kpa-1)*(order-cut+1)+1):(kpa*(order-cut+1)))=B0c;

            k(1,kpa)=kpa-1;
        end
        B=[BMr BMi];
    else
        % Se genera la matriz de Taylor-Fourier
        for kpa = 1:1:(2*(Nfrec)+1)
            for nI=1:1:(cut)
                for n = -Nh : 1: Nh
                    B0((n+Nh+1),nI)=double((n^(2*(nI-1)))*exp(j*n*w1*(kpa-1-Nfrec)));
                end
            end
            for nI=1:1:(order+1-cut)
                for n = -Nh : 1: Nh
                    B0c((n+Nh+1),nI)=double((n^(2*(nI-1)+1))*exp(j*n*w1*(kpa-1-Nfrec)));
                end    
            end
            BMr(1:N,((kpa-1)*cut+1):(kpa*cut))=B0;
            BMi(1:N,((kpa-1)*(order-cut+1)+1):(kpa*(order-cut+1)))=B0c;

            k(1,kpa)=kpa-1;
        end
        B=[BMr BMi];
    end

    Vent=ventana2d(Nh,alfa);

    est=(inv(B'*Vent*B)*B'*Vent);

    xx=linspace(-Nf/2,Nf/2, nmo );
   
    k1=k*cut+1;
    k2=((2*(Nfrec)+1)*cut)+(k*(order-cut+1))+1;
    k3=k*cut+2;
end

function [V] = ventana2d(Nh,alfa)
    N=2*Nh+1;
    M_Vent=zeros(1,(2*Nh)+1);		
    V=zeros(2*Nh+1,2*Nh+1);	
    M_Vent=(kaiser(N,alfa))';      
    for nI=1:1:((2*Nh)+1)			
        M_Vent(nI,:)=double(kaiser(N,alfa))';     	    
        V(nI,nI)=M_Vent(nI,nI)^(1);   
    end
end



