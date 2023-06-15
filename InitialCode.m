%% average data
clear;clc
sigma = 0.1; %noise parameter

A = [1 2; -3 2]; %stoichiometry example matrix
mub = [1; 1];  % average RHS measurements


mux = A\mub;   %expected solution to Ax=b

Sigma1 = diag([sigma^2, sigma^2*10]);  %Covariance matrix of measurement noise

Sigma2 = pinv(A)*Sigma1*pinv(A)';  %resulting covariance on x

[a, b, c, d] = deal(Sigma1(1,1),Sigma1(1,2),Sigma1(2,1),Sigma1(2,2));
det = a*d-b*c;

[a_,b_,c_,d_] = deal(Sigma2(1,1),Sigma2(1,2),Sigma2(2,1),Sigma2(2,2));
det_ = a_*d_-b_*c_;


alpha = 0.05; %p-value level
cst = chi2inv(1-alpha,length(mub)); %chi-2 quantile at p=(1-alpha) and degrees of freedom k = #variables

intervalb = [0,2];%plotting interval for measurements b
intervalx = [-0.25,0.25,0,1]; %plotting intervals for solutions x



%% simulation
bhat = [];
xhat = [];

subplot(121);title('measurements b');
hold on; axis equal;axis off;
subplot(122);title('solutions x');
hold on; axis equal;axis off;

%simulate 1000 samples from measurements distribution; associated solutions

for i =1:1000
    bhat(:,end+1) = mvnrnd(mub,Sigma1);

    xhat(:,end+1) = A\bhat(:,end);

end


subplot(121);
scatter(bhat(1,:),bhat(2,:),'b.');%mark observed measurements
scatter(mub(1),mub(2),'rx'); %mark mean measurement

%plot "cst"-levelset of measurement Mahalanobis distance; shoudl contain
%1-alpha of the observed samples bhat
fimplicit(@(x,y) (d*(x-mub(1)).^2 - b*(x-mub(1)).*(y-mub(2)) - c*(y-...
mub(2)).*(x-mub(1)) + a*(y-mub(2)).^2)/det - cst, intervalb,'r');


subplot(122);
scatter(xhat(1,:),xhat(2,:),'b.'); %mark observed solutions
scatter(mux(1),mux(2),'rx'); %mark mean solutions

% plot "cst"-levelset of solution Mahalanobis distance; should contain
% 1-alpha of the observed solutions xhat
fimplicit(@(x,y) (d_*(x-mux(1)).^2 - b_*(x-mux(1)).*(y-mux(2)) - c_*(y-...
mux(2)).*(x-mux(1)) + a_*(y-mux(2)).^2)/det_ - cst, intervalx,'r');
