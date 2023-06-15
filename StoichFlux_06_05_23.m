clear;clc;

%reading in the stoichimetry matix and converting to a double  
stoich_o = readtable('10237_2018_1068_MOESM1_ESM.xlsx');
stoich = removevars(stoich_o,{'Var1'});

stoich = table2array(stoich);

%reading in the flux matrix converting to a double
flux_o = readtable('10237_2018_1068_MOESM2_ESM.xlsx');
flux = removevars(flux_o,{'Var1'});

flux_15_0 = table2array(flux(:,1));
flux_15_15 = table2array(flux(:,2));

%compute the accumulations from the fluxes
c_15_0 = stoich*flux_15_0;
c_15_15 = stoich*flux_15_15;


%covariance of b - or the covariance matrix of metabolic intensity
%measurements
sigma=ones(length(c_15_0),1)*10^2;
e_b = diag(ones(length(c_15_0),1)*10^2).^2;



%% Compute the singular value decomposition of Stoich to determine the reaction rates
[U,S,V] = svd(stoich);


S_1 = S;
for i=1:size(S,2)
    %update the diagonals
    S_1(i,i) = 1/S(i,i);
end
S_1=S_1';

%calculating sigma_x, and u_x
mux = pinv(S)*c_15_0;
e_x = V*S_1*U'*e_b*U*S_1'*V';


