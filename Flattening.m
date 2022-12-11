clear all, close all, clc                   % clear all old variables
%% Flattening
% load petrological and model data
load Z_numeric_data.mat     % X coordinate measured W
load W_numeric.mat          % Y coordinate measured W
load Z_numeric_data_H.mat   % X coordinate measured H
load H_numeric.mat          % Y coordinate measured H
load Z_numeric.mat          % x cooridnate model W & H
load W.mat                  % y coordinate model W
load H.mat                  % y coordinate model H
% paraemter
X                       = Z_numeric;
X_area                  = 2;                                                % change the flattening limit (x-coordinate)
Ind                     = find(X<X_area); 
strain                  = -0.1;                                             % = 10% - change to find best fit
X_flat                  = X;
X_flat(Ind)             = (1 + strain) * X(Ind);
X_off                   = X(Ind(end)) - (1 + strain) * X(Ind(end));                 
X_flat(Ind(end)+1:end)  = X(Ind(end)+1:end) - X_off;
% plot ====================================================================
figure(1)
plot(Z_numeric_data,W_numeric,'-b'), hold on                                % data W
plot(Z_numeric,W,'k-'), hold on                                             % model W
plot(Z_numeric_data_H,H_numeric,'-r'),hold on                               % data H
plot(Z_numeric,H,'-m'),hold on                                              % model H
plot(X_flat,W,'-g'), hold on
plot(X_flat,H,'--g')
plot([1 1]*X_area,[-0.1 1],'--k')
axis([0 5 -0.05 1.05])
xlabel('Z')                                                                 % normalized distance
ylabel('W & H')
legend('W data','W model','H data','H model','10% flat W','10% flat H','flat area')








