clear;clc;close all;

Tsw = 2;
N = [1000, 1000, 1000, 1000, 1000];
ws = [1, 2.25, 4];

U1 = gbngen(N(1), Tsw);
w1 = ws(1) * ones(N(1), 1);
Y1 = zeros(N(1), 1);
for i = 1:N(1)
    [dnum,dden] = Get_dnumden(w1(i));
    Y1(i+1,1) = dnum(2)*U1(i) - dden(2) * Y1(i, 1);
end

U2 = gbngen(N(2), Tsw);
w2 = linspace(ws(1), ws(2), N(2))';
Y2 = zeros(N(2), 1);
% Y2(1, 1) = Y1(end, 1);
Y2(1, 1) = 0;
for i = 1:N(2)
    [dnum,dden] = Get_dnumden(w2(i));
    Y2(i+1, 1) = dnum(2)*U2(i) - dden(2) * Y2(i,1);
end

U3 = gbngen(N(3), Tsw);
w3 = ws(2)*ones(N(3), 1);
Y3 = zeros(N(3), 1);
% Y3(1, 1) = Y2(end, 1);
for i = 1:N(3)
    [dnum,dden] = Get_dnumden(w3(i));
    Y3(i+1, 1) = dnum(2)*U3(i) - dden(2) * Y3(i,1);
end

U4 = gbngen(N(4), Tsw);
w4 = linspace(ws(2), ws(3), N(4))';
Y4 = zeros(N(4), 1);
Y4(1, 1) = Y3(end, 1);
for i = 1:N(4)
    [dnum,dden] = Get_dnumden(w4(i));
    Y4(i+1, 1) = dnum(2)*U4(i) - dden(2) * Y4(i,1);
end

U5 = gbngen(N(5), Tsw);
w5 = ws(3) * ones(N(5), 1);
Y5 = zeros(N(5), 1);
Y5(1, 1) = Y4(end, 1);
for i = 1:N(5)
    [dnum,dden] = Get_dnumden(w5(i));
    Y5(i+1, 1) = dnum(2)*U5(i) - dden(2) * Y5(i,1);
end

% process Y
Y1 = Y1(1:N(1), 1);
Y2 = Y2(1:N(2), 1);
Y3 = Y3(1:N(3), 1);
Y4 = Y4(1:N(4), 1);
Y5 = Y5(1:N(5), 1);
% add noise 
% ratio = 0.03;
% Y_1 = Add_noise(Y_1, ratio);
% Y_2 = Add_noise(Y_2, ratio);
% Y_3 = Add_noise(Y_3, ratio);
% Y_4 = Add_noise(Y_4, ratio);
% Y_5 = Add_noise(Y_5, ratio);


da_Num = sum(N);
U = [U1;U2;U3;U4;U5];
W = [w1;w2;w3;w4;w5];
Y = [Y1;Y2;Y3;Y4;Y5];
% Y = zeros(da_Num, 1);
% calculate system output for a input u 
% for i = 1:da_Num
%     [dnum,dden] = Get_dnumden(W(i));
%     Y(i+1,1) = dnum(2)*U(i) - dden(2) * Y(i,1);
% end
% Y = Y(1:da_Num,1);

% PY = var(Y);
% n = randn(da_Num, 1);
% v = filter([1], [1, -0.9], n);
% Pv = var(v);
% Pn_desired = 0.03 * PY;
% v = v * sqrt(Pn_desired/Pv);
% SNR = 10 * log10(PY/Pn_desired);
% Y = Y + v;

% THoe01 = oe([Y(2001:2500,1),U_1],[1,1,1]); 
%% save data 
save('../FitModel/data/u.mat', 'U');
save('../FitModel/data/w.mat', 'W');
save('../FitModel/data/y.mat', 'Y');

%% load data 
clear;clc;close all;
U = load('../FitModel/data/u.mat').U;
W = load('../FitModel/data/w.mat').W;
Y = load('../FitModel/data/Y.mat').Y;
[da_Num, ~] = size(Y);
%% model Identification
phi_1 = [];
phi_2 = [];
phi_3 = [];
k = [1,1.5,1.9,2.25,2.7,3.4,4];
for i = 1:da_Num
    phi_1 = [phi_1;[1,W(i),abs(W(i)-k(2))^3,abs(W(i)-k(3))^3,abs(W(i)-k(4))^3,abs(W(i)-k(5))^3,abs(W(i)-k(6))^3]];
end
phi_2 = phi_1;
phi_3 = phi_1;

y1_hat =  Get_yhat(1,U);
y2_hat =  Get_yhat(2.25,U);
y3_hat =  Get_yhat(4,U);

PHI_1 = phi_1.*y1_hat;
PHI_2 = phi_2.*y2_hat;
PHI_3 = phi_3.*y3_hat;
PHI = [PHI_1,PHI_2,PHI_3];
beta = ((PHI'*PHI)\PHI')*Y;
w_step = (1:0.01:4)';
alpha_1 = beta(1) + beta(2)*w_step + beta(3)*abs(w_step-k(2)).^3 + beta(4)*abs(w_step-k(3)).^3 +...
    + beta(5)*abs(w_step-k(4)).^3 + beta(6)*abs(w_step-k(5)).^3 + beta(7)*abs(w_step-k(6)).^3;

alpha_2 = beta(1+7) + beta(2+7)*w_step + beta(3+7)*abs(w_step-k(2)).^3 + beta(4+7)*abs(w_step-k(3)).^3 +...
    + beta(5+7)*abs(w_step-k(4)).^3 + beta(6+7)*abs(w_step-k(5)).^3 + beta(7+7)*abs(w_step-k(6)).^3;

alpha_3 = beta(1+14) + beta(2+14)*w_step + beta(3+14)*abs(w_step-k(2)).^3 + beta(4+14)*abs(w_step-k(3)).^3 +...
    + beta(5+14)*abs(w_step-k(4)).^3 + beta(6+14)*abs(w_step-k(5)).^3 + beta(7+14)*abs(w_step-k(6)).^3;

figure(1);
plot(w_step,alpha_1,'-b',w_step,alpha_2,'--k',w_step,alpha_3,'-.r');
legend('Alpha_1','Alpha_2','Alpha_3')

Ynn = load("./data/nn.mat").y';
figure(2);
Nsim = 150;

work_pt = 1.5;
idx = find(w_step == work_pt);

True_STP = Get_stepresponse(work_pt,Nsim);
work_1 = Get_stepresponse(1,Nsim);
work_2 = Get_stepresponse(2.25,Nsim);
work_3 = Get_stepresponse(4,Nsim);
[w_x_1l,w_x_2l,w_x_3l] = Trivial_Interpolation(work_pt);
Tri_STP = w_x_1l*work_1 + w_x_2l*work_2 + w_x_3l*work_3;
LPV_STP = alpha_1(idx,1)*work_1 + alpha_2(idx,1)*work_2 + alpha_3(idx,1)*work_3;
plot(0:Nsim,[0;True_STP],'-b',0:Nsim,[0;Tri_STP],'--k',0:Nsim,[0;LPV_STP],'-.r', 0:Nsim, Ynn, '-*g');
legend('True','Linear','LPV', 'LSTM')

figure(3) 
error_Tri_STP = sum((True_STP - Tri_STP).^2, 'all');
error_LPV = sum((True_STP - LPV_STP).^2, 'all');
error_LSTM = sum(((True_STP - Ynn(2:end)).^2), 'all');
e = [error_Tri_STP, error_LPV, error_LSTM];
bar(e);
set(gca, 'XTickLabel', {'Trivial', 'LPV', 'LSTM'});
text(1:length(e), e, num2str(e'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');