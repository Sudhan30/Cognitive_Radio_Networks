clc;
clear all;
Pf = 0:0.01:1;
% Pf = Pf.^2;

SNR_dB= -12; % SNR in dB
gam=power(10,SNR_dB/10); % SNR ratio
N=1000; % No of samples
f=6e6; % carrier freq for modn
bw = 6e6; % bandwidth of interest
noise_figure = 11; % noise figure in dB (as specified in IEEE 802.22)
noise_power = -174+10*log(bw)+noise_figure; %noise power in dBm conversion
noise_pow=double((1e-3)*power(10,noise_power/10));% received power in watts

[bpsk time]=bpskmod(f); % Function call bpskmod
sig_pow = noise_pow*gam; %required signal power for the specified SNR
sig_dbw=10*log10(sig_pow); % signal power in dB
amp=double(sqrt(2*sig_pow)); % signal amp from signal power
sig=amp*bpsk; % bpsk signal with the required amplitude
% noise_pow = 1e-6;
% % Pd calculation using formula
%   sig(2001:3000)=0;


for i = 1:length(Pf) 
 lamda(i) = noise_pow*(N+(sqrt(2*N))*qfuncinv(Pf(i)));
    v1 = lamda(i)-(N*noise_pow*(1+gam))
    v2 = sqrt(2*N*(noise_pow^2)*((1+gam)^2))
    Pd(i) = qfunc(v1/v2);
end

% run=1000;
% for i = 1:length(Pf)
%     count = 0;
%     for j = 1:run        
%         sig_nse=awgn(sig(1:N),SNR_dB,sig_dbw);
%         noisesim=sum(sig_nse.^2)
%         
%         if (noisesim>lamda(i))
%             count=count+1;
%         end
%     end
%     
%     pdsim(i)=count/run;
% end

plot(Pf,Pd,'m','linewidth',1.5);axis ([0 1 0 1]);grid on;hold on;
% plot(Pf,pdsim,'r','linewidth',1.5);axis ([0 1 0 1]);
set(0,'DefaultAxesFontSize', 13);
set(gca,'FontSize',10,'fontWeight','bold');

xlabel('Probability of FALSE ALARM');
ylabel('Probability of DETECTION');
title('PROBABILITY OF FALSE ALARM VS DETECTION');
legend('Theoretical','Simulation','Location','North');

% Pd
% pdsim