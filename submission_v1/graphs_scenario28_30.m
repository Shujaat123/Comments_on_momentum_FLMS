clc
clear all
clf
loc=0;
label_ind=10;

for i=28:30%1:27

loc=loc+1;
load(strcat('Results/comp_LMS_mFLMS',int2str(i),'.mat'))

if(loc==1)
    plot(10*log10(mNWD_LMS(:)),'-','linewidth',2)
end
hold on 
plot(10*log10(mNWD_mFLMS(:)),'--','linewidth',2)
title(strcat('Noise level \sigma^2 = ',sprintf(' %.2f',(noise_level.^2))))
xlabel('No. of iterations')
ylabel('Fitness "\delta" (dB)')
grid minor
set(gca,'FontSize',13)


ylim([-25 0])
xlim([0 size(mNWD_mFLMS,2)])

if(loc==1)
leg_text{1}=['LMS_{\eta=',num2str(eta_LMS),'}'];
end
leg_text{loc+1}=['mFLMS_{f=',num2str(f_mFLMS),',\alpha=',num2str(alpha_mFLMS),',\eta=\eta_f=',num2str(eta_mFLMS),'}'];

if(loc==3)
    loc=0;
    h=legend(leg_text,'location','northeast');
    saveas(gcf,strcat('Figures\MSD_neu_vs_interations',int2str(label_ind),'.png'),'png')
%     saveas(gcf,strcat('Figures\MSD_neu_vs_interations',int2str(label_ind),'.eps'),'psc2')
    clf
    label_ind=label_ind+1;
end

end

