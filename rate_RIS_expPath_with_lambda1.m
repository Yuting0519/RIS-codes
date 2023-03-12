clear all
clc

%% Parameter Setup
tic
MC = 10^4; 
itrCh = 10^4; 
r0 = 0; %radius of inner circle
r1 = 10;   %radius of environment circle 
Cen_Xpos = 0; 
Cen_Ypos = 0;  
% ExpNum_Txs = 100; %number of TXs  
% dens = (ExpNum_Txs)/(pi*r1^2-pi*r0^2); 
lambda = 0.1:0.1:2; 
d = 1.2; 
S_Xpos = -d; S_Ypos = 0; 
D_Xpos = d; D_Ypos = 0;
N = 32; 
shape = N*pi^2/(16-pi^2); 
scale = (16-pi^2)/4/pi; 
pdB=0; 
p=10.^(0.1*pdB); 
% n0dB = -64; n0=10^(0.1*n0dB); 
n0=1;
alpha = 1.037; 
% thrdB = 5; thr = 10^(0.1*thrdB); 
% thr = 10^5; 
tp_opt = zeros(MC,length(lambda)); 
% Pout_prod = zeros(MC,length(pdB)); 
tp_min = zeros(MC,length(lambda)); 
% Pout_rad = zeros(MC,length(pdB)); 
tp_minMax = zeros(MC,length(lambda)); 
tp_minO = zeros(MC,length(lambda));
prob_analy = zeros(1,length(lambda));  
% Simulation
%every time step, one molecule is relased 
for mc = 1:MC
    for i = 1:length(lambda)
        dens=lambda(i);
        ExpNum_Txs = dens*(pi*r1^2-pi*r0^2); 
        PPNum_Txs = poissrnd(ExpNum_Txs);
        % generate HPPP positions of TXs and intial positions of molecules
        theta = rand(PPNum_Txs,1)*(2*pi);
        r = sqrt((r1^2-r0^2)*rand(PPNum_Txs,1)+r0^2);
        R_Xpos = Cen_Xpos + r.*cos(theta); %RX is at the center of circle
        R_Ypos = Cen_Ypos + r.*sin(theta);
        S_R = sqrt((R_Xpos-S_Xpos).^2+(R_Ypos-S_Ypos).^2);
        R_D = sqrt((R_Xpos-D_Xpos).^2+(R_Ypos-D_Ypos).^2);
        % distances from RIS to origin
        O_R = sqrt((R_Xpos-Cen_Xpos).^2+(R_Ypos-Cen_Ypos).^2);
        % choose the IR has minimum sum of SR+RD
        sum_SRD = S_R+R_D;
        minSum_SRD=min(sum_SRD);
        % choose the IR has mininal SR, RD
        [~, ind_min]=min(min(S_R,R_D));
        sum2_SRD=S_R(ind_min)+R_D(ind_min);
        % choose min-max
        [~, ind_minMax]=min(max(S_R,R_D));
        sum3_SRD=S_R(ind_minMax)+R_D(ind_minMax);
        % choose closet to O
        [~, ind_minO]=min(O_R);
        sum5_SRD=S_R(ind_minO)+R_D(ind_minO);

        %A = gamrnd(shape,scale,1,itrCh); % A generated from gamma
        % A generated from product Rayleigh    
        j=sqrt(-1);   
        h =  abs((1/sqrt(2))*(randn(N,itrCh)+ j*randn(N,itrCh)));  
        g =  abs((1/sqrt(2))*(randn(N,itrCh)+ j*randn(N,itrCh)));
        A = sum(h.*g);  
        tp_opt(mc,i) = mean(log2(1+p*A.^2/exp(alpha*minSum_SRD)/n0));
        tp_min(mc,i) = mean(log2(1+p*A.^2/exp(alpha*sum2_SRD)/n0));
        tp_minMax(mc,i) = mean(log2(1+p*A.^2/exp(alpha*sum3_SRD)/n0));
        tp_minO(mc,i) = mean(log2(1+p*A.^2/exp(alpha*sum5_SRD)/n0));
    end
end
tp1 = sum(tp_opt,1)/MC; 
tp3 = sum(tp_min,1)/MC;
tp4 = sum(tp_minMax,1)/MC;
tp6 = sum(tp_minO,1)/MC;

% y = log(p.*scale^2*gamma(2+shape)/thr/n0/gamma(shape))/alpha; 
% 
% %prob_analy = exp(-2*dens*(sqrt(1/d^4)*(d^4*ellipticE(asin(1),(y.^2/d^4)) ...
%     %+(y.^2-d^4).*ellipticF(asin(1),(y.^2/d^4)))));
% for i = 1:length(y)
%     if y(i) <= 2*d 
%         prob_analy(i)=1; 
%     else
% f = 1/8*(-1).^floor((pi-2.*angle(y(i))+angle(-4*d^2+y(i)^2))/2/pi)*pi.*y(i).*sqrt(-4*d^2+y(i)^2); 
% prob_analy(i) = exp(-2*dens*f);
%     end 
% end 
% tp_exp_PL=[lambda',tp1',tp6',tp3',tp4']
% save(sprintf('TPvsLam_exp_PL_P%d_N%d.mat',pdB,N),'tp_exp_PL');


figure = figure('PaperSize',[20.5 29],'Color',[1 1 1]);
axes = axes('Parent',figure,'FontSize',14,'FontName','Times New Roman');
plot(lambda, tp1, 'ro','Linewidth',1.2,'markers',7)
hold on
plot(lambda, tp3, 'bo','Linewidth',1.2,'markers',7)
plot(lambda, tp4, 'go','Linewidth',1.2,'markers',7)
plot(lambda, tp6, 'ko','Linewidth',1.2,'markers',7)
% plot(dens, prob_analy, 'r-','Linewidth',1.2,'markers',7)
xlabel('\lambda','FontSize',15,'FontName','Times New Roman');
ylabel('Average Rate','FontSize',15,'FontName','Times New Roman');
set(gca,'Fontsize',15);
grid on
box on;
set(gcf,'color','w');
set(gca,'FontName','Times New Roman');
% legend('optimal','mininal prod','miniMin','minMax','random','analytical');
legend('Optimum','Min-Min','Min-Max','Mid-Point','Analytical, Optimum');



