clear
clc
rng('shuffle')%set up the random seed
u_v=[200];% two possible values of change-points
K_d=1;% two clusters
m=40; % # of drivers
l1=0.25;% intensity rate before the change-point
l2=0.1; % intensity rate after the change-point
[z,Nj,C,tau_true_index]=latent_simu_f(u_v,m,l1,l2,K_d);%___________data simulation end_______________
% You can also make plots, I will show you next time
filename1=['censortm' num2str(1) '.csv'];
filename2=['eventtm' num2str(1) '.csv'];
csvwrite(filename1,C);
dlmwrite(filename2,z,'precision',10)