% clc;
clear all;
Pf = 0.1;
SNR_dB= -20:1:0; % SNR in dB
N=1000; 
noise_pow = 1;

m = 0.1:0.01:4;

% average mean and variance for the past L samples
M = 5;
L = 5;

%ROC curves using analytical expressions

for j = 1:length(SNR_dB) 
    gam=power(10,SNR_dB(j)/10); % SNR ratio
    %mean and variance of the pth summing detector
    for i = 1:length(m)
        p = m(i);
        mu_avg = ((M/L)*(1+gam)*noise_pow)+ (((L-M)/L)*noise_pow);
        sig_avg = ((M/(L^2))* (2/N)*(((1+gam)^2)*(noise_pow^2))+(((L-M)/(L^2))*(2/N)*(noise_pow^2)));
        mu0 = (2^(p/2))*gamma((p+1)/2)/sqrt(pi);
        var0 = (2^p)*(gamma((2*p+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi)))/sqrt(pi);
        mu1 = (2^(p/2))*((1+gam)^(p/2))*gamma((p+1)/2)/sqrt(pi);
        var1 = (2^p)*((1+gam)^(p))*(gamma((2*p+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi)))/sqrt(pi);
        lamda = (qfuncinv(Pf)*(sqrt(var0)/sqrt(N)))+mu0;
        Pd(i) = qfunc(sqrt(N)*(lamda-mu1)/sqrt(var1));
        Pfaied(i) = Pf+(Pf*((1-Pf)*qfunc(((lamda)-mu_avg)/sqrt(sig_avg))));
        Pdied(i) = Pd(i)+(Pd(i)*((1-Pd(i))*qfunc(((lamda)-mu_avg)/sqrt(sig_avg))));
    end
    Pe = Pfaied+(1-Pdied);
    [val pos] = min(Pe);
    opt_m(j)=m(pos);
    Pe_min(j)=val;
%     [val pos] = max(Pdied);
%     opt_m(j)=m(pos);
%     Pd_max(j)=val;
end

plot(SNR_dB,Pe_min,'-om','linewidth',1.5);
hold on;

for i = 1:length(SNR_dB) 
       gam=power(10,SNR_dB(i)/10); % SNR ratio
    %mean and variance of the pth summing detector
        p = 2;
        mu_avg = ((M/L)*(1+gam)*noise_pow)+ (((L-M)/L)*noise_pow);
        sig_avg = ((M/(L^2))* (2/N)*(((1+gam)^2)*(noise_pow^2))+(((L-M)/(L^2))*(2/N)*(noise_pow^2)));
        mu0 = (2^(p/2))*gamma((p+1)/2)/sqrt(pi);
        var0 = (2^p)*(gamma((2*p+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi)))/sqrt(pi);
        mu1 = (2^(p/2))*((1+gam)^(p/2))*gamma((p+1)/2)/sqrt(pi);
        var1 = (2^p)*((1+gam)^(p))*(gamma((2*p+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi)))/sqrt(pi);
        lamda = (qfuncinv(Pf)*(sqrt(var0)/sqrt(N)))+mu0;
        Pd1(i) = qfunc(sqrt(N)*(lamda-mu1)/sqrt(var1));
        Pfaied1(i) = Pf+(Pf*((1-Pf)*qfunc(((lamda)-mu_avg)/sqrt(sig_avg))));
        Pdied1(i) = Pd1(i)+(Pd1(i)*((1-Pd1(i))*qfunc(((lamda)-mu_avg)/sqrt(sig_avg))));
        Pe1 = Pfaied1+(1-Pdied1);
        
%        [val pos] = max(Pdied);
%     opt_m(j)=m(pos);
%     Pd_max(j)=val;
end
plot(SNR_dB,Pe1,'-*r')
hold on
plot(SNR_dB,Pec,'--b');
xlabel('SNR (dB)');
ylabel('Probability of error');
legend('proposed','IED','TED')