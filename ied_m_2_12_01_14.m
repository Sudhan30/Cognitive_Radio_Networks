clc;
clear all;
Pf = 0:0.01:1;
% Pf = Pf.^2;

SNR_dB= -10; % SNR in dB
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
% sig(1:2000)=0;

M =8;
L = 10;
mu_avg = ((M/L)*N*(1+gam)*noise_pow)+ (((L-M)/L)*N*noise_pow);
sig_avg = ((M/(L^2))* (2*N)*(((1+gam)^2)*(noise_pow^2))+(((L-M)/(L^2))*(2*N)*(noise_pow^2)));

for i = 1:length(Pf) 
    lamda(i) = (((qfuncinv(Pf(i)))*(sqrt(2*N)))+N)*noise_pow;
%   lamda(i) = noise_pow*(1+(sqrt(2/N))*qfuncinv(Pf(i)));
    v1 = sqrt((2*N*(gam+1)^2*noise_pow^2));
    v2 = lamda(i) - (N*(1+gam)*noise_pow);
    Pd(i) = qfunc(v2/v1);
%     Pfamed(i) = Pf(i)+(((1-Pf(i))*qfunc((lamda(i)-mu_avg)/sqrt(sig_avg))));
%     Pdmed(i) = Pd(i)+(((1-Pd(i))*qfunc((lamda(i)-mu_avg)/sqrt(sig_avg))));
    Pdied(i) = Pd(i)+(Pd(i)*((1-Pd(i))*qfunc((lamda(i)-mu_avg)/sqrt(sig_avg))));
     Pfied(i) = Pf(i)+(Pf(i)*((1-Pf(i))*qfunc((lamda(i)-mu_avg)/sqrt(sig_avg))));
%     Pdmed(i) = Pd(i)+(Pd(i)*((1-Pd(i))));
end

% run=1000;
% 
% for kk=1:L
%         sig_nse=awgn(sig(1:N),SNR_dB,sig_dbw);
%         noisesim=sum(sig_nse.^2);
%         ti_avg(kk)= noisesim;
% end
% 
% % ti_avg
% 
% for i = 1:length(Pf)
%     count =0;
% %     ti_avg(1:L)=0;
%     for j = 1:run       
%         mul=rem(j,L+1);
%         sig_nse=awgn(sig((N*mul)+1:(N*mul+N)),SNR_dB,sig_dbw);
%         noisesim=sum(sig_nse.^2);
%         
%         for jj=1:(L-1)
%             ti_avg(jj)=ti_avg(jj+1);
%         end
%         
%         ti_avg(L)=noisesim;
%         
%         if (noisesim>lamda(i))
%             count=count+1;
%         else
%             if (mean(ti_avg)>lamda(i))
%                 if (ti_avg(L-1)>lamda(i))
%                     count=count+1;
%             end
%         end
%         end
%     end
%  
%     pdsim(i)=count/run;
% end


plot(Pfied,Pdied,'k+-','linewidth',1.5);axis ([0 1 0 1]);grid on;hold on;
plot(Pf,Pd)
% plot(Pf,pdsim,'+-m','linewidth',1.5);axis ([0 1 0 1]);
xlabel('Probability of FALSE ALARM');
ylabel('Probability of DETECTION');
title('PROBABILITY OF FALSE ALARM VS DETECTION');
