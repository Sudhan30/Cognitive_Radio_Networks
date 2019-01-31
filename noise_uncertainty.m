clc;
clear all;
Pf = 0:0.1:1;
% Pf = Pf.^2;

SNR_dB= -10; % SNR in dB
gam=power(10,SNR_dB/10); % SNR ratio
N=1000; % No of samples
f=6e6; % carrier freq for modn
bw = 6e6; % bandwidth of interest
noise_figure = 11; % noise figure in dB (as specified in IEEE 802.22)
noise_power = -174+10*log(bw)+noise_figure; %noise power in dBm conversion
noise_pow=double((1e-3)*power(10,noise_power/10));% received power in watts
% noise_pow = 1;
[bpsk time]=bpskmod(f); % Function call bpskmod
sig_pow = noise_pow*gam; %required signal power for the specified SNR
sig_dbw=10*log10(sig_pow); % signal power in dB
amp=double(sqrt(2*sig_pow)); % signal amp from signal power
sig=amp*bpsk; % bpsk signal with the required amplitude
% noise_pow = 1e-6;
% % Pd calculation using formula
%   sig(2001:3000)=0;
del=0;
rho=10^(del/10);
    for i = 1:length(Pf) 
        for theor=1:1000
     a=noise_pow*(1/rho);
     b=noise_pow*rho;
     lamda_low(i) =a *(N+(sqrt(2*N))*qfuncinv(Pf(i)));
     lamda_high(i) =b*(N+(sqrt(2*N))*qfuncinv(Pf(i)));  
     r = a + (b-a).*rand(1,1);
     v1 = lamda_high(i)-(N*r*(1+gam));
     v2 = sqrt(2*N*(r^2)*((1+gam)^2));
     Pd(theor) = qfunc(v1/v2);
        end
    pd(i)=mean(Pd);
    end
    
% run=1000;
%  aa=SNR_dB-del;
%  ba=SNR_dB+del;
% for i = 1:length(Pf)
%     count = 0;
%     for j = 1:run        
%       
%         noise_un= aa+ (ba-aa).*rand(1,1);
%         sig_nse=awgn(sig(1:N),noise_un,sig_dbw);
%         noisesim=sum(sig_nse.^2);
%         
%         if (noisesim>lamda_high(i))
%             count=count+1;
%         end
%     end
%     
%     pdsim(i)=count/run;
% end

plot(Pf,pd,'r','linewidth',1.5);axis ([0 1 0 1]);grid on;hold on;
% plot(Pf,pdsim,'r','linewidth',1.5);axis ([0 1 0 1]);
set(0,'DefaultAxesFontSize', 13);
set(gca,'FontSize',10,'fontWeight','bold');

xlabel('Probability of FALSE ALARM');
ylabel('Probability of DETECTION');
title('PROBABILITY OF FALSE ALARM VS DETECTION');
% legend('Theoretical','Simulation','Location','North');
% 
% Pd
% pdsim