clear all, close all, clc           % clear old variables
%% Visualizing
load Data_PLOT_W_nl
load Data_PLOT_H_nl
load Data_PLOT_SZ
%% ------------------------------------
% plot for W
figure(2) 
sym         = {'-k','-b','-r','c','g'};
sym2        = {'--k','--b','--r','--c','--g'};

% plot W
for pp=1:length(Gamma_PLOT_W_nl(:,1))
    subplot(411)
    plot(Z/(SZ_gamma(end)/2),Gamma_PLOT_W_nl(pp,:),sym{pp})
    hold on
end
for pp=1:length(Gamma_PLOT_W_nl(:,1)) % to plot the shear zone thickness
    plot([1 1]*SZ_gamma(pp)/SZ_gamma(end),[0 2],sym2{pp})
    hold on
end
legend('\gamma = 0','\gamma = 0.25','\gamma = 0.5','\gamma = 0.75','\gamma = 1')
xlabel('X/k')
ylabel('W')
xlim([0 5])
ylim([0 1])
title('inflow of H_2O - eclogitization')

% plot normalized displacment
for pp=1:length(Dis_PLOT_W_nl(:,1))
    subplot(412)
    plot(Z/(SZ_gamma(end)/2), 2*(Dis_PLOT_W_nl(pp,:)-mean(Dis_PLOT_W_nl(pp,:))) ,sym{pp})
    hold on
end
for pp=1:length(Gamma_PLOT_W_nl(:,1)) % to plot the shear zone thickness
    plot([1 1]*SZ_gamma(pp)/SZ_gamma(end),[0 2],sym2{pp})
    hold on
end
xlabel('X/k')
ylabel('U')  
xlim([0 5])
ylim([0 1.1])
title('Displacement')

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
title('Viscosity with W')

% % plot variable diffusivity DW       % comment line 57 - 66 if DW is constant
% for pp=1:length(DW_var_PLOT_nl(:,1))
%     subplot(414)
%     plot(Z,DW_var_PLOT_nl(pp,:),sym{pp})
%     hold on
% end
% xlabel('Z')
% ylabel('DH')
% xlim([0 100])
% title('Diffusivity of W')

%% ------------------------------------------------------------------------
% plot for H
figure(3) 
sym         = {'-k','-b','-r','c','g'};
sym2        = {'--k','--b','--r','--c','--g'};
% plot H
for pp=1:length(Gamma_PLOT_H_nl(:,1))
    subplot(311)
    plot(Z,Gamma_PLOT_H_nl(pp,:),sym{pp})
    hold on
end
for pp=1:length(Gamma_PLOT_H_nl(:,1))
    plot([1 1]*SZ_gamma(pp)/SZ_gamma(end),[0 2],sym2{pp})
    hold on
end
legend('\gamma = 0','\gamma = 0.25','\gamma = 0.5','\gamma = 0.75','\gamma = 1')
xlabel('X/k')
ylabel('H')
xlim([0 5])
ylim([0 1])
title('influx of H - hydration of NAMs')

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
title('Viscosity with H')

% % plot variable diffusivity DH       % comment line 102 - 111 if DH is constant
% for pp=1:length(DH_var_PLOT_nl(:,1))
%     subplot(313)
%     plot(Z,DH_var_PLOT_nl(pp,:),sym{pp})
%     hold on
% end
% xlabel('Z')
% ylabel('DH')
% xlim([0 100])
% title('Diffusivity of H')

