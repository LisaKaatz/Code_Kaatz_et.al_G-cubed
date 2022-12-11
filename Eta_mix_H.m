function [A_bulk_H]=Eta_mix(Xr,A1,A3)

nb_mat = 2;   % number of considered materials
%-----------------------------------------------------
% INPUT INPUT INPUT
%-----------------------------------------------------
A = zeros(1,nb_mat);
f = zeros(1,nb_mat);

%====================================================
% % Material 1 (Matrix): Strong
A(1) = A1;
f(1) = 1.0-Xr;
% Material "3" : Transition
A(2) = A3;
f(2) = Xr;
a    = [2 2];
n    = [1 1];

Prod1 = 1;
sum_down =0;
sum_n = 0;
Prod2 = 1;

for i=1:nb_mat
    sum_down = sum_down + f(i)*a(i);
end

for i=1:nb_mat
    Prod1 = Prod1*A(i)^(f(i)*a(i)/sum_down);
    sum_n = sum_n + f(i)*n(i)/(n(i)+1);
    Prod2 = Prod2*(n(i)/(n(i)+1))^(f(i)*a(i)*n(i)/sum_down);
end

A_bulk_H = Prod1*sum_n^(-1)*Prod2;

end
