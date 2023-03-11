clear variables, close all, clc              % clear variables etc.
%% General model geometry and resolution
L            = 400;                          % model height / sample length
nz           = 301;                          % number of grid points - use 601 for accurate, but slow simulation 
dz           = L/(nz-1);                     % grid spacing
Z            = -L/2:dz:L/2;                  % coordinate vector
%% Parameters for aqueous fluid inflow
DW           = 20;                           % diffusivity of aqueous fluid
W0           = 1;                            % initial concentration = 1 = 100% eclogite phases
W            = W0*zeros(1,nz);               % concentration array, for outside of the shear zone
I_W          = [(nz+1)/2-1:(nz+1)/2+1];      % initial sz width, area where W & H = 1
W(I_W)       = 1;                            % assign W, for inside of the shear zone
W_0          = W;                            % initial W array
Res_W        = ones(1,length(W));            % residual for W
D_R          = 4;                            % diffusivity ratio (DH/DW)
% % if DW is variable:
%    a_W          = 0.01;                      % factor to change diffusivity between 0 and 1
%    b_W          = 0.5;                       % factor to change DW based on W
%    W_ratio      = ((1+a_W)./(W+a_W)).^b_W;   % non-linear equation for W
%    DW_var       = W_ratio*DW;                % variable DW
%% Parameters for hydrogen diffusion 
DH           = D_R*DW;                       % diffusivity of hydrogen 
H0           = 1;                            % initial concentration = 1 = maximum water content
H            = H0*zeros(1,nz);               % concentration array, for outside of the shear zone
I_H          = I_W;                          
H(I_H)       = 1;                            % assign H, for inside of the shear zone
H_0          = H;                            % initial H array
Res_H        = ones(1,length(H));            % Residual for H
%% if DH is variable:
%     a_H          = 0.01;                     % factor to change diffusivity between 0 and 1
%     b_H          = 0.5;                      % factor to change DH based on H
%     H_ratio      = ((1+a_H)./(W+a_H)).^b_H;  % factor to change relation of non-linearity
%     DH_var       = H_ratio*DH;               % variable DH
%% Parameters for viscous shear deformation
Vx           = zeros(1,nz+1);                % initial velocity
Vx(1)        = -0.5; Vx(end) = 0.5;          % top and bottom boundary velocity
%% Parameters for viscosity 1 - strong
eta1         = 1;                            % strong viscosity
A1           = 1/(2*eta1);                   % preexponential factor for strong viscosity
%% Parameters for viscosity 2 - weak
eta2         = 1e-4;                         % weak viscosity
A2           = 1/(2*eta2);                   % preexponential factor for weak viscosity
%% Parameters for viscosity 3 - transition
eta_R        = 100;                          % viscosity ratio transition/weak 
eta3         = eta_R*eta2;                   % viscosity of transition
A3           = 1/(2*eta3);                   % preexponential factor for transition
eta_Huet     = ones(1,length(W));            % Initialization for viscosity array
%% Numerical parameters
I            = 2:length(Z)-1;                           % indices
dt           = dz^2/max([DW DH])/2;                     % time step                 
dt_pt        = dz^2/max([DW DH])/1e1;                   % iteration step            
Time         = [0];                                     % Time array
time         = 0;                                       % initial time
tolerance    = 3e-6;                                    % error for iterative solver - use 1e-7 for accurate, but slow simulation 
ti           = 0;                                       % initial time step counter  
%% Parameters to save
Displacement_tot = 0; Gamma_tot_s1     = 0;  Gamma_tot_s2     = 0;                                   
Dis_prof_tot     = 0; Strain_prof      = 0;  gamma_exit       = 0;
%     DW_var_nl   = 0;        % if DW is variable
%     DH_var_nl   = 0;        % if DH is variable
%% Initializing matrices for later visualization
save_plot_count  = 1;
% gamma
Gamma_PLOT_W_nl           = zeros(5,nz);
Gamma_PLOT_W_nl(1,:)      = W;
Gamma_PLOT_H_nl           = zeros(5,nz);
Gamma_PLOT_H_nl(1,:)      = H;
gamma_save_ref            = 0.25;           % to plot gamma within 5 steps 0 - 1
gamma_save                = 0.25;
% displacement
Dis_PLOT_W_nl             = zeros(5,nz);
Dis_PLOT_W_nl(1,:)        = zeros(1,nz);
Dis_PLOT_H_nl             = zeros(5,nz);
Dis_PLOT_H_nl(1,:)        = zeros(1,nz);
% strain
Strain_PLOT_W_nl          = zeros(5,nz);
Strain_PLOT_W_nl(1,:)     = zeros(1,nz);
Strain_PLOT_H_nl          = zeros(5,nz);
Strain_PLOT_H_nl(1,:)     = zeros(1,nz);
% % if DW or DH are variable         
%     DW_var_PLOT_nl            = zeros(5,nz);
%     DW_var_PLOT_nl(1,:)       = zeros(1,nz);
%     DH_var_PLOT_nl            = zeros(5,nz);
%     DH_var_PLOT_nl(1,:)       = zeros(1,nz);
% effective viscosity
Visc_PLOT_W_nl            = zeros(5,nz);
Visc_PLOT_H_nl            = zeros(5,nz);
%SZthickness
SZ_gamma = 0;
for ii=1:length(H)
                [A_bulk_H]        = Eta_mix_H(H(ii),A1,A3);
                eta_Huet_H(ii)    = 0.5*A_bulk_H^(-1);
end
A3_H = 1./(2*eta_Huet_H);
for ii=1:length(W)
                [A_bulk]        = Eta_mix_W(W(ii),A3_H(ii),A2);
                eta_Huet(ii)    = 0.5*A_bulk^(-1);
end
Visc_PLOT_H_nl(1,:)       = eta_Huet_H;
Visc_PLOT_W_nl(1,:)       = eta_Huet;
%% Start of time loop
while gamma_exit<1
    ti           = ti+1;                     % time step counter
    W_old        = W;                        % update W
    H_old        = H;                        % update H
    error        = 1e6;
    it           = 0;                        % initial interation counter
    while error>tolerance | it<30
        it       = it+1;
        % ===== Diffusion ======
        % linear diffusion of W
        Res_W(I)        = -(W(I)-W_old(I))/dt + ( DW*(W(I+1)-2*W(I)+W(I-1))/dz^2 );                           % Residual of W equation
        W(I)            = W(I) + dt_pt*Res_W(I);                                                              % iterative update W
        W(I_W)          = 1;                                                                                  % keep initial W in center
            % % if DW is variable 
            % W_ratio      = ((1+a_W)./(W+a_W)).^b_W;   % non-linearity equation for W
            % DW_var       = W_ratio*DW;                % variable DW
            % Flux         = -(DW_var(1:end-1) + DW_var(2:end))/2 .* (W(2:end) - W(1:end-1))/dz; 
            % Res_W(I)     = -(W(I)-W_old(I))/dt - ( Flux(2:end) - Flux(1:end-1))/dz;                            % Residual of W equation
            % W(I)         = W(I) + dt_pt*Res_W(I);                                                              % iterative update W
            % W(I_W)       = 1;                                                                                  % keep initial W in center
       % linear diffusion of H
        Res_H(I)        = -(H(I)-H_old(I))/dt + ( DH*(H(I+1)-2*H(I)+H(I-1))/dz^2 );                           % Residual of H equation
        H(I)            = H(I) + dt_pt*Res_H(I);                                                              % iterative update H
        H(I_H)          = 1;                                                                                  % keep initial H in center
%         % if DH is variable
%         H_ratio      = ((1+a_H)./(H+a_H)).^b_H;   % non-linearity equation for W
%         DH_var       = H_ratio*DH;                % variable DW
%         Flux_H          = -(DH_var(1:end-1) + DH_var(2:end))/2 .* (H(2:end) - H(1:end-1))/dz; 
%         Res_H(I)        = -(H(I)-H_old(I))/dt - ( Flux_H(2:end) - Flux_H(1:end-1))/dz;                        % Residual of H equation
%         H(I)            = H(I) + dt_pt*Res_H(I);                                                              % iterative update C
%         H(I_H)          = 1;                                                                                  % keep initial H in center                                                        
        % ===== Viscous flow ==== Viscosity mixing
        if mod(it,25)==1 | it==1
            for ii=1:length(H)
                [A_bulk_H]        = Eta_mix_H(H(ii),A1,A3);
                eta_Huet_H(ii)    = 0.5*A_bulk_H^(-1);
            end
        end
        % Effective viscosity of granulite is now an array as function of H 
        A3_H = 1./(2*eta_Huet_H);                               
        if mod(it,25)==1 | it==1
            for ii=1:length(W)
                [A_bulk]        = Eta_mix_W(W(ii),A3_H(ii),A2);
                eta_Huet(ii)    = 0.5*A_bulk^(-1);
            end
        end
        Eta         = eta_Huet;
        tau         = 2*Eta.*diff(Vx)/dz;                                   % shear stress
        res_Vx      = diff(tau)/dz;                                         % residual of force balance equation
        dt_visc     = dz^2/max(Eta(:))/4;                                   % iteration parameter
        Vx(2:end-1) = Vx(2:end-1) + dt_visc*res_Vx;                         % velocity update
        Res_W(I_W)  = 0;                                                    % in case central W is constant
        Res_H(I_H)  = 0;                                                    % in case central H is constant
        error_W     = max(abs(Res_W(I)));                                   % error for W
        error_H     = max(abs(Res_H(I)));                                   % error for H
        error_Vx    = max(abs(res_Vx));                                     % error for Vx
        error       = max([error_W error_H error_Vx]);                      % total error
    end
    time          = time+dt;                                                % time update
    Time          = [Time time];                                            % store time in time array
    if ti==1; Vx0 = (Vx(1:end-1)+Vx(2:end))/2; end                          % safe initial velocity
        Vxc       = (Vx(1:end-1)+Vx(2:end))/2;                              
        Strain_rate         = diff(Vx)/dz;                                  
        Strain              = Strain_rate*dt;
        Strain_prof         = Strain_prof + Strain;                         
        
        % Calculation of discplacement
        Dis_prof_inc        = cumsum(Strain*dz);                            
        Dis_prof_tot        = Dis_prof_tot + Dis_prof_inc;
        
        % SZ thickness with strain (> meanStrain)
        I_SZ                = find(Strain_prof > mean(Strain_prof));                                
        SZthickness(ti)     = Z(I_SZ(end))-Z(I_SZ(1));                    
        Gamma(ti)           = (Dis_prof_tot(I_SZ(end)) - mean(Dis_prof_tot)) / (SZthickness(ti)/2); 
        gamma_exit          = Gamma(end);
      
   % Output after timesteps
   if mod(ti,50)==1                                                         % save/show every 25th time step
        figure(1)
        subplot(241), hold off
        plot(Z,Vxc,'-k'), hold on, plot(Z,Vx0,'--k'), xlabel('X'), ylabel('Velocity')
        legend('Velocity','Initial velocity','location','northwest'), 
        title(['Error: ',num2str(error,3),' Iter.: ',num2str(it)])      
        
        % plot W & H
        subplot(242), hold off
        plot(Z,W,'-k'), hold on
        plot(Z,W_0,'--k'), hold on
        plot(Z,H,'-r')
        xlabel('X'), ylabel('W & H'), axis([-L/2 L/2 0 1]), legend('W','initial W','H bulk')
        title(['Time step: ',num2str(ti)]), hold on
        
        % plot effective viscosity
        subplot(243), hold off
        plot(Z,Eta,'-k'), hold on                                                   
        plot(Z,eta_Huet_H,'-r')                                                     
        xlabel('X'), ylabel('Viscosity'), axis([-L/2 L/2 0 1])
        set(gca,'yscale','log')
        grid on
        legend('\eta mix total','\eta mix strong')
                       
        % plot shear zone thickness
        subplot(245)
        plot(Time(2:end),SZthickness,'-r')
        legend('SZT (meanStrain)','location','southeast')
        xlabel('Time'), ylabel('Shear zone thickness')
        
        % gamma plotten
        subplot(246)
        plot(Time(2:end),Gamma,'b-')
        xlabel('Time'), ylabel('Gamma')
        
        % displacement plotten
        subplot(247), hold off
        plot(Z,Dis_prof_tot/(SZthickness(end)),'k-')
        ylabel('Normalized displacement'),xlabel('X')
        
        subplot(248)
        plot(Z,Strain_prof,'-k'), hold on
        plot([Z(1) Z(end)],[mean(Strain_prof) mean(Strain_prof)],'-r'), hold off
        legend('Strain','Mean strain')
        ylabel('Strain'),xlabel('X')

        set(gcf,'position',[405.8000 264.2000 1.0744e+03 420.0000])
        drawnow
   end   
   %% Storing into matrices for later visualization
    if Gamma(end)>gamma_save
        save_plot_count                     = save_plot_count + 1;
        Gamma_PLOT_W_nl(save_plot_count,:)  = W;
        Gamma_PLOT_H_nl(save_plot_count,:)  = H;
        Dis_PLOT_W_nl(save_plot_count,:)    = Dis_prof_tot/SZthickness(end);
        Visc_PLOT_W_nl(save_plot_count,:)   = Eta;
        Visc_PLOT_H_nl(save_plot_count,:)   = eta_Huet_H;
        Strain_PLOT_W_nl(save_plot_count,:) = Strain_prof;
        SZ_gamma(save_plot_count)           = SZthickness(end);
%             DW_var_PLOT_nl(save_plot_count,:)       = DW_var;               % if DW is variable
%             DH_var_PLOT_nl(save_plot_count,:)       = DH_var;               % if DH is variable
        gamma_save                          = gamma_save_ref * save_plot_count;
    end 
end
%% Saving the matrices a .mat files for visualization
save Data_PLOT_W_nl Gamma_PLOT_W_nl Dis_PLOT_W_nl Visc_PLOT_W_nl DW_var_PLOT_nl Z;           % add: DW_var_PLOT_nl if DW is variable
save Data_PLOT_H_nl Gamma_PLOT_H_nl Dis_PLOT_H_nl Visc_PLOT_H_nl Z;           % add: DH_var_PLOT_nl if DH is variable
save Data_PLOT_SZ   SZ_gamma Z; 

%% load/plot petrological data
load Cx.mat; load Cy.mat;                                                   % normalized eclogite-facies mineral assemblage
load bulkHZ.mat; load HPhases.mat;                                          % normalized water-content data

% INTERP on grid points ===================================================
Z_numeric           = Z/(SZthickness(end)/2);
Ind_data            = find(Cx(1) < Z_numeric & Z_numeric < Cx(end));        % interpolate model W-Data
Z_numeric_data      = Z_numeric(Ind_data);
W_numeric           = interp1(Cx,Cy,Z_numeric_data);
Ind_data_H          = find(HPhases(1) < Z_numeric & Z_numeric < HPhases(end));
Z_numeric_data_H    = Z_numeric(Ind_data_H);
H_numeric           = interp1(HPhases,BulkHZ,Z_numeric_data_H);
% save for later mechanical pure shear correction =========================
save Z_numeric.mat; save W_numeric.mat; save H_numeric.mat; save H.mat; save W.mat; save Z_numeric_data.mat; save Z_numeric_data_H.mat
%% Plot results at gamma = 1
figure(2)
subplot(211)
plot(Cx,Cy,'k-'); hold on                       % W data fit
plot(Z_numeric_data,W_numeric,'sk')             % INTERP data C
plot(HPhases,BulkHZ,'r-'); hold on              % H data fit 
plot(Z_numeric_data_H,H_numeric,'sr')           % INTERP data H 
plot(Z_numeric,W,'-*k')                         % W model 
plot(Z_numeric,H,'-*r')                         % H model 

% ERROR calculation =======================================================
Error_W_vec         = abs(W_numeric - W(Ind_data))./W_numeric;
Error_H_vec         = abs(H_numeric - H(Ind_data_H))./H_numeric;
Error_W_av          = mean(Error_W_vec); maxError_W_av = max(Error_W_vec);  minError_W_av = min(Error_W_vec);
Error_H_av          = mean(Error_H_vec); maxError_H_av = max(Error_H_vec);  minError_H_av = min(Error_H_vec);
Error_tot           = (Error_W_av + Error_H_av)/2;
Error_W             = Error_W_av*100;    maxError_W    = maxError_W_av*100; minError_W    = minError_W_av*100;
Error_H             = Error_H_av*100;    maxError_H    = maxError_H_av*100; minError_H    = minError_H_av*100;
Error_T             = Error_tot*100;
% ERROR ===================================================================

axis([0 L/2 0 1.0])
xlabel('X / (SZT/2)'),ylabel('W & H')                                             % normalized x-axis
legend('fit W','interp W','fit H bulk','interp H','W model','H model')
title(['results \gamma = 1','- Mean Error W: ',num2str(Error_W),'%', '- Max Error W: ',num2str(maxError_W),'%', '- Min Error W: ',num2str(minError_W),'%',...
                            '- Mean Error H: ',num2str(Error_H),'%', '- Max Error H: ',num2str(maxError_H),'%', '- Min Error H: ',num2str(minError_H),'%'])
xlim([0 5])
grid on

subplot(212)
plot(Z/(SZthickness(end)/2),Strain_prof,'-xk'), hold on
plot(Z/(SZthickness(end)/2),(Dis_prof_tot-mean(Dis_prof_tot))/(SZthickness(end)/2),'-*r'), hold on          
hei         = max([max(Strain_prof) max((Dis_prof_tot-mean(Dis_prof_tot))/(SZthickness(end)/2))]); 
plot([1 1],[0 hei*1.1],'--k')
xlim([0 5]);ylim([0 hei*1.1])
title('Strain and displacement profile')
xlabel('X / (SZT/2)'), ylabel('Strain & displacement')
legend('Strain','Displacement / (SZT/2)','Shear zone boundary')
grid on
set(gcf,'Position',[81 67.4000 1.1864e+03 516])