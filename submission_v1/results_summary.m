clc
clear all
close all

Results = [];
for i=19:27
    load(strcat('Results/comp_LMS_mFLMS',int2str(i),'.mat'))
    
    Results=[Results; (noise_level.^2), alpha_mFLMS, f_mFLMS, eta_mFLMS, theta_mFLMS, mse(theta,theta_mFLMS)*1e10; ...
        (noise_level.^2), alpha_mFLMS, f_mFLMS, eta_LMS,theta_LMS, mse(theta,theta_LMS)*1e10];
end

Results=[Results;0, 0, 0, 0,theta, 0];
Results = round(10000*Results)/10000;
Results(:,13) = Results(:,13)/1e10;


