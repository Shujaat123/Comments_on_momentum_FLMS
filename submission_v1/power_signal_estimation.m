clc; clear all; 
for scenario=1:27
    close all;

%% Power Signal parameter estimation using mFLMS and LMS
runs=100;
%% Change Noise Level to 0.3
%% Change fractional power to 0.25
if (scenario==1) % Tested
    N=1:8000;   
    eta_LMS = 2.9e-3;
    noise_level = sqrt(0.3);
    f_mFLMS = 0.25;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;

elseif (scenario==2)% Tested
    N=1:8000;
    eta_LMS = 4.4e-3;
    noise_level = sqrt(0.3);
    f_mFLMS = 0.25;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
           
elseif (scenario==3)% Tested
    N=1:8000;
    eta_LMS = 1.15e-2;
    noise_level = sqrt(0.3);
    f_mFLMS = 0.25;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;
%% Change fractional power to 0.50
elseif (scenario==4) % Tested
    N=1:8000;   
    eta_LMS = 2.9e-3;
    noise_level = sqrt(0.3);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;
    
elseif (scenario==5) % Tested
    N=1:8000;
    eta_LMS = 4.4e-3;
    noise_level = sqrt(0.3);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
    
elseif (scenario==6) % Tested
    N=1:8000;
    eta_LMS = 1.15e-2;
    noise_level = sqrt(0.3);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;
    
%% Change fractional power to 0.75
elseif (scenario==7) % Tested
    N=1:8000;   
    eta_LMS = 2.9e-3;
    noise_level = sqrt(0.3);
    f_mFLMS = 0.75;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;
    
elseif (scenario==8) % Tested
    N=1:8000;   
    eta_LMS = 4.4e-3;
    noise_level = sqrt(0.3);
    f_mFLMS = 0.75;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
    
elseif (scenario==9) % Tested
    N=1:8000;   
    eta_LMS = 1.15e-2;
    noise_level = sqrt(0.3);
    f_mFLMS = 0.75;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;
    
%% Change Noise Level to 0.6
%% Change fractional power to 0.25
elseif (scenario==10) % Tested
    N=1:8000;   
    eta_LMS = 2.8e-3;
    noise_level = sqrt(0.6);
    f_mFLMS = 0.25;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;
    
elseif (scenario==11) % Tested
    N=1:8000;   
    eta_LMS = 4.4e-3;
    noise_level = sqrt(0.6);
    f_mFLMS = 0.25;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
    
elseif (scenario==12) % Tested
    N=1:8000;   
    eta_LMS = 1.15e-2;
    noise_level = sqrt(0.6);
    f_mFLMS = 0.25;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;
    
%% Change fractional power to 0.50
elseif (scenario==13)% Tested
    N=1:8000;   
    eta_LMS = 2.8e-3;
    noise_level = sqrt(0.6);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;
    
elseif (scenario==14)% Tested
    N=1:8000;   
    eta_LMS = 4.4e-3;
    noise_level = sqrt(0.6);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
    
elseif (scenario==15)
    N=1:8000;   
    eta_LMS = 1.15e-2;
    noise_level = sqrt(0.6);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;
    
%% Change fractional power to 0.75
elseif (scenario==16)% Tested
    N=1:8000;   
    eta_LMS = 2.8e-3;
    noise_level = sqrt(0.6);
    f_mFLMS = 0.75;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;
    
elseif (scenario==17) % Tested
    N=1:8000;   
    eta_LMS = 4.4e-3;
    noise_level = sqrt(0.6);
    f_mFLMS = 0.75;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
    
elseif (scenario==18)% Tested
    N=1:8000;   
    eta_LMS = 1.15e-2;
    noise_level = sqrt(0.6);
    f_mFLMS = 0.75;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;
    
%% Change Noise Level to 0.9
%% Change fractional power to 0.25
elseif (scenario==19)% Tested
    N=1:8000;   
    eta_LMS = 2.9e-3;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.25;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;
    
elseif (scenario==20)% Tested
    N=1:8000;   
    eta_LMS = 4.4e-3;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.25;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
    
elseif (scenario==21)% Tested
    N=1:8000;   
    eta_LMS = 1.15e-2;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.25;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;
    
%% Change fractional power to 0.50
elseif (scenario==22)% Tested
    N=1:8000;   
    eta_LMS = 2.9e-3;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;
    
elseif (scenario==23)% Tested
    N=1:8000;   
    eta_LMS = 4.4e-3;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
    
elseif (scenario==24)% Tested
    N=1:8000;   
    eta_LMS = 1.15e-2;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;
    
%% Change fractional power to 0.75
elseif (scenario==25)% Tested
    N=1:8000;   
    eta_LMS = 2.9e-3;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.75;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;
    
elseif (scenario==26)% Tested
    N=1:8000;   
    eta_LMS = 4.4e-3;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.75;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
    
elseif (scenario==27)% Tested
    N=1:8000;
    eta_LMS = 1.15e-2;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.75;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;

elseif (scenario==28)% Tested
    N=1:16000;
    eta_LMS = 1e-3;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.2;
elseif (scenario==29)% Tested
    N=1:16000;
    eta_LMS = 1e-3;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.50;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.5;
elseif (scenario==30)% Tested
    N=1:16000;
    eta_LMS = 1e-3;
    noise_level = sqrt(0.9);
    f_mFLMS = 0.5;
    eta_mFLMS = 1e-3;
    eta_f_mFLMS = eta_mFLMS/gamma(2-f_mFLMS);
    alpha_mFLMS = 0.8;
end


%%---------------------------------------%%
theta = [1.8 2.9 4 2.5 0.95 0.8 0.76 1.1];
omega = [0.07 0.5 2 1.6];

mtheta_recov_LMS = 0;
mNWD_LMS = 0;
mMSE_LMS = 0;

mtheta_recov_mFLMS=0; % mean value of recovered theta initialized to be zero
mNWD_mFLMS=0;         % mean normalized weight deviation initialized to be zero
mMSE_mFLMS=0;         % mean mean squared error initialized to be zero

%% Input signal
for n = 1:length(N)
X_input(:,n) = ([sin(omega(1)*n) cos(omega(1)*n) sin(omega(2)*n) cos(omega(2)*n) sin(omega(3)*n) cos(omega(3)*n) sin(omega(4)*n) cos(omega(4)*n)]);
end

%% Ideal desired signal
for n=1:length(N)
    r(n) = theta(1)*sin(omega(1)*n +theta(5)) + theta(2)*sin(omega(2)*n + theta(6))+ theta(3)*sin(omega(3)*n + theta(7)) + theta(4) *sin(omega(4)*n +theta(8));
end 

for k=1:runs
    d = r + noise_level*randn(size(r)); % noisy desired signal
    %% parameters for mFLMS
    U_mFLMS = zeros(1,length(theta));
    W_mFLMS = randn(1,length(theta));
    e_mFLMS = 0;
    temp_mFLMS = 0;
    NWD_mFLMS = zeros(1,length(N));
    
    %% parameters for LMS 
    U_LMS = zeros(1,length(theta));
    W_LMS = W_mFLMS;
    e_LMS = 0;
    NWD_LMS = zeros(1,length(N));
    
for n=1:length(N)
    %% mLMS Algorithm
    U_mFLMS=X_input(:,n)';   
    Y_mFLMS = W_mFLMS*U_mFLMS';
    e_mFLMS = d(n) - Y_mFLMS;
    temp_mFLMS = alpha_mFLMS*temp_mFLMS + eta_mFLMS*e_mFLMS*U_mFLMS + eta_f_mFLMS*e_mFLMS*U_mFLMS.*(abs(W_mFLMS).^(1-f_mFLMS));
    W_mFLMS = W_mFLMS + temp_mFLMS;
    
theta_recov_mFLMS(1)=norm(W_mFLMS(1:2));
theta_recov_mFLMS(5)=atan(W_mFLMS(2)/W_mFLMS(1));
theta_recov_mFLMS(2)=norm(W_mFLMS(3:4));
theta_recov_mFLMS(6)=atan(W_mFLMS(4)/W_mFLMS(3));
theta_recov_mFLMS(3)=norm(W_mFLMS(5:6));
theta_recov_mFLMS(7)=atan(W_mFLMS(6)/W_mFLMS(5));
theta_recov_mFLMS(4)=norm(W_mFLMS(7:8));
theta_recov_mFLMS(8)=atan(W_mFLMS(8)/W_mFLMS(7));


%% LMS Algorithm
    
    U_LMS=X_input(:,n)';   
    Y_LMS = W_LMS*U_LMS';
    e_LMS = d(n) - Y_LMS;
    W_LMS = W_LMS + eta_LMS* e_LMS*U_LMS;
    
theta_recov_LMS(1)=norm(W_LMS(1:2));
theta_recov_LMS(5)=atan(W_LMS(2)/W_LMS(1));
theta_recov_LMS(2)=norm(W_LMS(3:4));
theta_recov_LMS(6)=atan(W_LMS(4)/W_LMS(3));
theta_recov_LMS(3)=norm(W_LMS(5:6));
theta_recov_LMS(7)=atan(W_LMS(6)/W_LMS(5));
theta_recov_LMS(4)=norm(W_LMS(7:8));
theta_recov_LMS(8)=atan(W_LMS(8)/W_LMS(7));

NWD_mFLMS(n) = norm(theta-theta_recov_mFLMS)./norm(theta);
NWD_LMS(n) = norm(theta-theta_recov_LMS)./norm(theta);

end
mNWD_LMS = mNWD_LMS + NWD_LMS;
mtheta_recov_LMS=mtheta_recov_LMS+theta_recov_LMS;

mNWD_mFLMS = mNWD_mFLMS + NWD_mFLMS;
mtheta_recov_mFLMS=mtheta_recov_mFLMS+theta_recov_mFLMS;

end

mNWD_LMS=mNWD_LMS/runs;
mNWD_mFLMS=mNWD_mFLMS/runs;
theta_mFLMS = mtheta_recov_mFLMS/runs;
theta_LMS = mtheta_recov_LMS/runs;

plot(10*log10(mNWD_mFLMS),'k','linewidth',2)
hold on 
plot(10*log10(mNWD_LMS),'r','linewidth',2)
title('Momentum FLMS/LMS')
xlabel('No. of iterations')
ylabel('Fitness "\delta" (dB)')
h=legend('mFLMS','LMS');
grid minor
set(h,'FontSize',13)
set(gca,'FontSize',13)
saveas(gcf,strcat('MSD_neu_vs_interations',int2str(scenario),'.png'),'png')


[theta_mFLMS;theta_LMS;theta]

save(strcat('Results/comp_LMS_mFLMS',int2str(scenario),'.mat'));

end
