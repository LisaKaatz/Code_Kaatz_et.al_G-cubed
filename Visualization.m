clear all, close all, clc           % clear old variables
%% Visualizing
load Data_PLOT_W_nl
load Data_PLOT_H_nl

% plot for W
figure(2) 
sym         = {'-k','-b','-r','c','g'};
% plot W
for pp=1:length(Gamma_PLOT_W_nl(:,1))
    subplot(411)
    plot(Z,Gamma_PLOT_W_nl(pp,:),sym{pp})
    hold on
end
legend('\gamma = 0','\gamma = 0.25','\gamma = 0.5','\gamma = 0.75','\gamma = 1')
xlabel('Z')
ylabel('W')
xlim([0 100])
ylim([0 1])
title('propagation of W')
% plot normalized displacment
for pp=1:length(Dis_PLOT_W_nl(:,1))
    subplot(412)
    plot(Z,Dis_PLOT_W_nl(pp,:),sym{pp})
    hold on
end
xlabel('Z')
ylabel('U')  
xlim([0 100])
ylim([0 1])
title('norm. displacement with W')
% plot effective viscosity
for pp=1:length(Visc_PLOT_W_nl(:,1))
    subplot(413)
    plot(Z,Visc_PLOT_W_nl(pp,:),sym{pp})
    hold on
    set(gca,'yscale','log')
end
xlabel('Z')
ylabel('\eta eff')
xlim([0 100])
title('eff viscosity with W')
% plot variable diffusivity DW       % comment line 43-52 if DW is constant
for pp=1:length(DW_var_PLOT_nl(:,1))
    subplot(313)
    plot(Z,DW_var_PLOT_nl(pp,:),sym{pp})
    hold on
end
xlabel('Z')
ylabel('DH')
xlim([0 100])
title('diffusivity with H')

% ------------------------------------------------------------------------
% plot for H
figure(3) 
sym         = {'-k','-b','-r','c','g'};
% plot H
for pp=1:length(Gamma_PLOT_H_nl(:,1))
    subplot(311)
    plot(Z,Gamma_PLOT_H_nl(pp,:),sym{pp})
    hold on
end
legend('\gamma = 0','\gamma = 0.25','\gamma = 0.5','\gamma = 0.75','\gamma = 1')
xlabel('Z')
ylabel('H')
xlim([0 100])
ylim([0 1])
title('propagation of H')
% plot effective viscosity
for pp=1:length(Visc_PLOT_H_nl(:,1))
    subplot(312)
    plot(Z,Visc_PLOT_H_nl(pp,:),sym{pp})
    hold on
    set(gca,'yscale','log')
end
xlabel('Z')
ylabel('\eta eff')
xlim([0 100])
title('eff viscosity with H')
% plot variable diffusivity DH       % comment line 71-72 if DH is constant
for pp=1:length(DH_var_PLOT_nl(:,1))
    subplot(313)
    plot(Z,DH_var_PLOT_nl(pp,:),sym{pp})
    hold on
end
xlabel('Z')
ylabel('DH')
xlim([0 100])
title('diffusivity with H')
