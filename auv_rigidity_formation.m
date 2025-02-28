clc;
clear;
close all;
addpath('./utilis/','./header/')
set(0,'defaultfigurecolor','w')
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');
%======================================================================%
%======================================================================%

err_model=0;
m=116;
Iz=13.1;
X_udot=-167.6*(1+err_model);
Y_vdot=-477.2*(1+err_model);
N_rdot=-15.9*(1+err_model);
Xu=26.9*(1+err_model);
Yv=35.8*(1+err_model);
Nr=3.5*(1+err_model);
Du=241.3*(1+err_model);
Dv=503.8*(1+err_model);
Dr=76.9*(1+err_model);

Mx=m-X_udot;
My=m-Y_vdot;
Mpsi=Iz-N_rdot;
disturbance=[0;0;0];
%======================================================================%
%======================================================================%

coef = [m; Iz; X_udot; Y_vdot; N_rdot; Xu; Yv; Nr; Du; Dv; Dr];
ndof = 3;

x1 = 0;
y1 = 0;

x2 = -7.42;
y2 = 6.7;  

x3 = 4.08;
y3 = 9.13;

x4 = 9.94;
y4 = -1.06;

x5 = 2.07;
y5 = -9.78;

x6 = -8.67;
y6 = -4.99;

X1_0 = [x1       ;y1       ;  0 + pi  ; 0; 0; 0];
X2_0 = [x2 - 2.76;y2 + 6.14;  0 + pi/2; 0; 0; 0];
X3_0 = [x3 + 7.3 ;y3 - 1.4 ;  0 + pi/3; 0; 0; 0];
X4_0 = [x4 + 2.3 ;y4 - 4.2 ;  0 - pi/2; 0; 0; 0];
X5_0 = [x5 + 10.4;y5 - 6.6 ;  0 + pi/4; 0; 0; 0];
X6_0 = [x6 - 5.4 ;y6 - 10.2;  0 - pi/3; 0; 0; 0];


U_0 = [0;0;0];

auv1 = AUV(coef, ndof, X1_0, U_0);
auv2 = AUV(coef, ndof, X2_0, U_0);
auv3 = AUV(coef, ndof, X3_0, U_0);
auv4 = AUV(coef, ndof, X4_0, U_0);
auv5 = AUV(coef, ndof, X5_0, U_0);
auv6 = AUV(coef, ndof, X6_0, U_0);

internal_auv1 = AUV(coef, ndof, X1_0, U_0);
internal_auv2 = AUV(coef, ndof, X2_0, U_0);
internal_auv3 = AUV(coef, ndof, X3_0, U_0);
internal_auv4 = AUV(coef, ndof, X4_0, U_0);
internal_auv5 = AUV(coef, ndof, X5_0, U_0);
internal_auv6 = AUV(coef, ndof, X6_0, U_0);

%======================================================================%
%======================================================================%

Adj = [0,1,1,1,1,1;
       1,0,1,0,0,1;
       1,1,0,1,0,0;
       1,0,1,0,1,0;
       1,0,0,1,0,0;
       1,1,0,0,0,0];

ref = [0,    1,    1,    1,    1,    1;
       1,    0, 1.18,    0,    0, 1.18;
       1, 1.18,    0, 1.18,    0,    0;
       1,    0, 1.18,    0, 1.18,    0;
       1,    0,    0, 1.18,    0,    0;
       1, 1.18,    0,    0,    0,    0]*10;



traj = [0,0,0,0,0,0];
%==========================================================================
%==========================================================================

dt = 0.1;

Tstep = 150;

X1plus = X1_0;
X2plus = X2_0;
X3plus = X3_0; 
X4plus = X4_0; 
X5plus = X5_0; 
X6plus = X6_0;

%======================================================================%
%======================================================================%

nx = length(X1_0);
nu = length(U_0);
[Xa1, Xa2, Xa3, Xa4, Xa5, Xa6] = deal(zeros(nx,Tstep+1));
[Ua1, Ua2, Ua3, Ua4, Ua5, Ua6] = deal(zeros(nu,Tstep));
[Va1, Va2, Va3, Va4, Va5, Va6] = deal(zeros(nu,Tstep));
[Vap1, Vap2, Vap3, Vap4, Vap5, Vap6] = deal(zeros(nu,Tstep));

Xa1(:,1) = X1_0; 
Xa2(:,1) = X2_0; 
Xa3(:,1) = X3_0; 
Xa4(:,1) = X4_0;
Xa5(:,1) = X5_0;
Xa6(:,1) = X6_0;

k_bs.v = 1e-5;
k_bs.a = 1e3;
k_bs.psi_v = 1;
k_bs.psi_a = 1;

upperbound_bs = 1000;
lowerbound_bs = -1000;

bs1 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
bs2 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
bs3 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
bs4 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
bs5 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
bs6 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);

%======================================================================%
%======================================================================%

for i=1:1:Tstep
    state = [X1plus,X2plus,X3plus,X4plus,X5plus,X6plus];
    [tau1,Vdot1] = bs1.calc_control(ref,traj,1,state);
    [tau2,Vdot2] = bs2.calc_control(ref,traj,2,state);
    [tau3,Vdot3] = bs3.calc_control(ref,traj,3,state);
    [tau4,Vdot4] = bs4.calc_control(ref,traj,4,state);
    [tau5,Vdot5] = bs5.calc_control(ref,traj,5,state);
    [tau6,Vdot6] = bs6.calc_control(ref,traj,6,state);
    
    Ua1(:,i) = tau1;
    Ua2(:,i) = tau2;
    Ua3(:,i) = tau3;
    Ua4(:,i) = tau4;
    Ua5(:,i) = tau5;
    Ua6(:,i) = tau6;

    Va1(:,i) = Vdot1(1);
    Va2(:,i) = Vdot2(1);
    Va3(:,i) = Vdot3(1);
    Va4(:,i) = Vdot4(1);
    Va5(:,i) = Vdot5(1);
    Va6(:,i) = Vdot6(1);

    Vap1(:,i) = Vdot1(2);
    Vap2(:,i) = Vdot2(2);
    Vap3(:,i) = Vdot3(2);
    Vap4(:,i) = Vdot4(2);
    Vap5(:,i) = Vdot5(2);
    Vap6(:,i) = Vdot6(2);

    auv1.advance(tau1, disturbance, dt);
    X1plus = auv1.X;
    Xa1(:,i+1) = X1plus;

    auv2.advance(tau2, disturbance, dt);
    X2plus = auv2.X;
    Xa2(:,i+1) = X2plus;

    auv3.advance(tau3, disturbance, dt);
    X3plus = auv3.X;
    Xa3(:,i+1) = X3plus;

    auv4.advance(tau4, disturbance, dt);
    X4plus = auv4.X;
    Xa4(:,i+1) = X4plus;

    auv5.advance(tau5, disturbance, dt);
    X5plus = auv5.X;
    Xa5(:,i+1) = X5plus;

    auv6.advance(tau6, disturbance, dt);
    X6plus = auv6.X;
    Xa6(:,i+1) = X6plus;
end


%%
%======================================================================%
%======================================================================%
blue   = '#5F97D2';
lgreen = '#B1CE46';
red    = '#D76364';
yellow = '#F1D77E';
green  = '#63E398';
purple = '#7E2F8E';
pink   = '#A2142F';
cblue  = '#5F9EA0';
close all;
figure(1);

plot(Xa1(1,end), Xa1(2,end),'*','color',blue,'LineWidth',2);
hold on;
plot(Xa2(1,end), Xa2(2,end),'*','color',lgreen,'LineWidth',2);
hold on;
plot(Xa3(1,end), Xa3(2,end),'*','color',red,'LineWidth',2);
hold on;
plot(Xa4(1,end), Xa4(2,end),'*','color',yellow,'LineWidth',2);
hold on;
plot(Xa5(1,end), Xa5(2,end),'*','color',green,'LineWidth',2);
hold on;
plot(Xa6(1,end), Xa6(2,end),'*','color',purple,'LineWidth',2);
hold on;

b1 = plot([Xa1(1,1),Xa2(1,1)],[Xa1(2,1),Xa2(2,1)],'--','color',cblue);
hold on;
plot([Xa1(1,1),Xa3(1,1)],[Xa1(2,1),Xa3(2,1)],'--','color',cblue);
hold on;
plot([Xa1(1,1),Xa4(1,1)],[Xa1(2,1),Xa4(2,1)],'--','color',cblue);
hold on;
plot([Xa1(1,1),Xa5(1,1)],[Xa1(2,1),Xa5(2,1)],'--','color',cblue);
hold on;
plot([Xa1(1,1),Xa6(1,1)],[Xa1(2,1),Xa6(2,1)],'--','color',cblue);
hold on;
plot([Xa2(1,1),Xa3(1,1)],[Xa2(2,1),Xa3(2,1)],'--','color',cblue);
hold on;
plot([Xa3(1,1),Xa4(1,1)],[Xa3(2,1),Xa4(2,1)],'--','color',cblue);
hold on;
plot([Xa4(1,1),Xa5(1,1)],[Xa4(2,1),Xa5(2,1)],'--','color',cblue);
hold on;
plot([Xa5(1,1),Xa6(1,1)],[Xa5(2,1),Xa6(2,1)],'--','color',cblue);
hold on;
plot([Xa2(1,1),Xa6(1,1)],[Xa2(2,1),Xa6(2,1)],'--','color',cblue);
hold on;

b2 = plot([Xa1(1,end),Xa2(1,end)],[Xa1(2,end),Xa2(2,end)],'--','color',pink);
hold on;
plot([Xa1(1,end),Xa3(1,end)],[Xa1(2,end),Xa3(2,end)],'--','color',pink);
hold on;
plot([Xa1(1,end),Xa4(1,end)],[Xa1(2,end),Xa4(2,end)],'--','color',pink);
hold on;
plot([Xa1(1,end),Xa5(1,end)],[Xa1(2,end),Xa5(2,end)],'--','color',pink);
hold on;
plot([Xa1(1,end),Xa6(1,end)],[Xa1(2,end),Xa6(2,end)],'--','color',pink);
hold on;
plot([Xa2(1,end),Xa3(1,end)],[Xa2(2,end),Xa3(2,end)],'--','color',pink);
hold on;
plot([Xa3(1,end),Xa4(1,end)],[Xa3(2,end),Xa4(2,end)],'--','color',pink);
hold on;
plot([Xa4(1,end),Xa5(1,end)],[Xa4(2,end),Xa5(2,end)],'--','color',pink);
hold on;
plot([Xa5(1,end),Xa6(1,end)],[Xa5(2,end),Xa6(2,end)],'--','color',pink);
hold on;
plot([Xa2(1,end),Xa6(1,end)],[Xa2(2,end),Xa6(2,end)],'--','color',pink);
hold on;

a1=plot(Xa1(1,:), Xa1(2,:),'color',blue,'LineWidth',2);
hold on;
a2=plot(Xa2(1,:), Xa2(2,:),'color',lgreen,'LineWidth',2);
hold on;
a3=plot(Xa3(1,:), Xa3(2,:),'color',red,'LineWidth',2);
hold on;
a4=plot(Xa4(1,:), Xa4(2,:),'color',yellow,'LineWidth',2);
hold on;
a5=plot(Xa5(1,:), Xa5(2,:),'color',green,'LineWidth',2);
hold on;
a6=plot(Xa6(1,:), Xa6(2,:),'color',purple,'LineWidth',2);
hold on;
legend([a1,a2,a3,a4,a5,a6,b1,b2],...
{'{BSC UUV1}','{BSC UUV2}', '{BSC UUV3}', '{BSC UUV4}','{BSC UUV5}','{BSC UUV6}','{Initial formation}','{Final formation}'}, ...
'NumColumns',4,'location','northeast','Interpreter','latex');
legend('boxoff');

axis square
xlabel('\boldmath{$x[m]$}', 'Interpreter', 'latex','FontSize',14);
ylabel('\boldmath{$y[m]$}', 'Interpreter', 'latex','FontSize',14);
set(gca,'FontSize',12);
set(gcf,'unit','normalized','position', [0.2,0.2,0.6,0.6]);

% axis equal

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
t = tiledlayout(3,1);
nexttile;
a1=plot(0.1:0.1:Tstep/10,Ua1(1,:),'color',blue,'LineWidth',2);
hold on;
a2=plot(0.1:0.1:Tstep/10,Ua2(1,:),'color',lgreen,'LineWidth',2);
hold on;
a3=plot(0.1:0.1:Tstep/10,Ua3(1,:),'color',red,'LineWidth',2);
hold on;
a4=plot(0.1:0.1:Tstep/10,Ua4(1,:),'color',yellow,'LineWidth',2);
hold on;
a5=plot(0.1:0.1:Tstep/10,Ua5(1,:),'color',green,'LineWidth',2);
hold on;
a6=plot(0.1:0.1:Tstep/10,Ua6(1,:),'color',purple,'LineWidth',2);
legend([a1,a2,a3,a4,a5,a6],...
{'BSC UUV1','BSC UUV2', 'BSC UUV3', 'BSC UUV4','BSC UUV5','BSC UUV6'},'NumColumns',3,'location','northeast');
legend('boxoff');
set(gca,'FontSize',12);
xlabel('\bf{time}\boldmath{$[s]$}', 'Interpreter', 'latex','FontSize',14);
ylabel('\textbf{control input} \boldmath{$\tau^{x}[N]$}', 'Interpreter', 'latex','FontSize',14);
nexttile;

a1=plot(0.1:0.1:Tstep/10,Ua1(2,:),'color',blue,'LineWidth',2);
hold on;
a2=plot(0.1:0.1:Tstep/10,Ua2(2,:),'color',lgreen,'LineWidth',2);
hold on;
a3=plot(0.1:0.1:Tstep/10,Ua3(2,:),'color',red,'LineWidth',2);
hold on;
a4=plot(0.1:0.1:Tstep/10,Ua4(2,:),'color',yellow,'LineWidth',2);
hold on;
a5=plot(0.1:0.1:Tstep/10,Ua5(2,:),'color',green,'LineWidth',2);
hold on;
a6=plot(0.1:0.1:Tstep/10,Ua6(2,:),'color',purple,'LineWidth',2);
xlabel('\bf{time}\boldmath{$[s]$}', 'Interpreter', 'latex','FontSize',14);
ylabel('\textbf{control input} \boldmath{$\tau^{y}[N]$}', 'Interpreter', 'latex','FontSize',14);
set(gca,'FontSize',12);
nexttile;

a1=plot(0.1:0.1:Tstep/10,Ua1(3,:),'color',blue,'LineWidth',2);
hold on;
a2=plot(0.1:0.1:Tstep/10,Ua2(3,:),'color',lgreen,'LineWidth',2);
hold on;
a3=plot(0.1:0.1:Tstep/10,Ua3(3,:),'color',red,'LineWidth',2);
hold on;
a4=plot(0.1:0.1:Tstep/10,Ua4(3,:),'color',yellow,'LineWidth',2);
hold on;
a5=plot(0.1:0.1:Tstep/10,Ua5(3,:),'color',green,'LineWidth',2);
hold on;
a6=plot(0.1:0.1:Tstep/10,Ua6(3,:),'color',purple,'LineWidth',2);
xlabel('\bf{time}\boldmath{$[s]$}', 'Interpreter', 'latex','FontSize',14);
h = ylabel('\textbf{control input} \boldmath{$\tau^{r}[N\cdot m]$}', 'Interpreter', 'latex','FontSize',14);
set(gca,'FontSize',12);

set(gcf,'unit','normalized','position', [0.2,0.2,0.5,0.5]);

%%

auv1 = AUV(coef, ndof, X1_0, U_0);
auv2 = AUV(coef, ndof, X2_0, U_0);
auv3 = AUV(coef, ndof, X3_0, U_0);
auv4 = AUV(coef, ndof, X4_0, U_0);
auv5 = AUV(coef, ndof, X5_0, U_0);
auv6 = AUV(coef, ndof, X6_0, U_0);

%======================================================================%
%======================================================================%
N = 4;  
M = nx + nu;
[u1_0, u2_0, u3_0, u4_0, u5_0, u6_0] = deal(zeros(nu*N,1));


U1_max = 1e3;
U2_max = 1e3;
U3_max = 1e3;

U_max0 = [U1_max;U2_max;U3_max];
U_min0 = [-U1_max;-U2_max;-U3_max];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Li = diag([5e3 1e2 1e2 1e2]);
Qi = diag([0.6 1 1 1 1 1])*3e3;
Ri = diag([1 1 1])*1e-2;
Qfi = diag([0.6 1 1 1 1 1])*3e3;
weights = {Li,Qi,Ri,Qfi};
traj = zeros(9,4);
traj(5:6,:)=0;
auxiliary_controller1 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
auxiliary_controller2 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
auxiliary_controller3 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
auxiliary_controller4 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
auxiliary_controller5 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);
auxiliary_controller6 = rigidity_bs_controller(k_bs,coef,upperbound_bs,lowerbound_bs,Adj);

dlmpc1 =rigidity_dlmpc_controller(N,internal_auv1,auxiliary_controller1,weights,U_max0,U_min0,Adj);
dlmpc2 =rigidity_dlmpc_controller(N,internal_auv2,auxiliary_controller2,weights,U_max0,U_min0,Adj);
dlmpc3 =rigidity_dlmpc_controller(N,internal_auv3,auxiliary_controller3,weights,U_max0,U_min0,Adj);
dlmpc4 =rigidity_dlmpc_controller(N,internal_auv4,auxiliary_controller4,weights,U_max0,U_min0,Adj);
dlmpc5 =rigidity_dlmpc_controller(N,internal_auv5,auxiliary_controller5,weights,U_max0,U_min0,Adj);
dlmpc6 =rigidity_dlmpc_controller(N,internal_auv6,auxiliary_controller6,weights,U_max0,U_min0,Adj);

%======================================================================%
%======================================================================%
[X1, X2, X3, X4, X5, X6] = deal(zeros(nx,Tstep+1));
[U1, U2, U3, U4 ,U5, U6] = deal(zeros(nu,Tstep));
state_pred = zeros(nx,N,6);

for i = 1:6
    x_pred = eval(['X', num2str(i),'_0']);
    for j = 1:4
        state_pred(:, j, i) = x_pred;
    end
end

X1(:,1) = X1_0;
X2(:,1) = X2_0; 
X3(:,1) = X3_0;
X4(:,1) = X4_0;
X5(:,1) = X5_0;
X6(:,1) = X6_0;

t = 0;
state = [X1_0,X2_0,X3_0,X4_0,X5_0,X6_0];

tStart = tic;
for i=1:1:Tstep
    fprintf("T=%d\n",i);
    % AUV1 DLMPC
    tic % 计时 
    [u1, X1_pred] = dlmpc1.calc_control(ref,1,traj,state,state_pred,u1_0,dt);
    state_pred(:,:,1) = X1_pred;
    u1_actual = u1(1:nu,1);
    U1(:,i) = u1_actual;
    auv1.advance(u1_actual, disturbance, dt);
    
    X1(:,i+1) = auv1.X;
    state(:,1) = auv1.X;
    u1_0 = dlmpc1.calc_initial_guess(ref,traj,1,state,dt);
    toc

    % AUV2 DLMPC
    tic
    [u2, X2_pred] = dlmpc2.calc_control(ref,2,traj,state,state_pred,u2_0,dt);
    state_pred(:,:,2) = X2_pred;
    u2_actual = u2(1:nu,1);
    U2(:,i) = u2_actual;
    auv2.advance(u2_actual, disturbance, dt);
    X2(:,i+1) = auv2.X;
    state(:,2) = auv2.X;
    u2_0 = dlmpc2.calc_initial_guess(ref,traj,2,state,dt);
    toc

    % AUV3 DLMPC
    tic
    [u3, X3_pred] = dlmpc3.calc_control(ref,3,traj,state,state_pred,u3_0,dt);
    state_pred(:,:,3) = X3_pred;
    u3_actual = u3(1:nu,1);
    U3(:,i) = u3_actual;
    auv3.advance(u3_actual, disturbance, dt);
    X3(:,i+1) = auv3.X;
    state(:,3) = auv3.X;
    u3_0 = dlmpc3.calc_initial_guess(ref,traj,3,state,dt);
    toc

    % AUV4 DLMPC
    tic
    [u4, X4_pred] = dlmpc4.calc_control(ref,4,traj,state,state_pred,u4_0,dt);
    state_pred(:,:,4) = X4_pred;
    u4_actual = u4(1:nu,1);
    U4(:,i) = u4_actual;
    auv4.advance(u4_actual, disturbance, dt);
    X4(:,i+1) = auv4.X;
    state(:,4) = auv4.X;
    u4_0 = dlmpc4.calc_initial_guess(ref,traj,4,state,dt);  
    toc

    % AUV5 DLMPC
    tic
    [u5, X5_pred] = dlmpc5.calc_control(ref,5,traj,state,state_pred,u5_0,dt);
    state_pred(:,:,5) = X5_pred;
    u5_actual = u5(1:nu,1);
    U5(:,i) = u5_actual;
    auv5.advance(u5_actual, disturbance, dt);
    X5(:,i+1) = auv5.X;
    state(:,5) = auv5.X;
    u5_0 = dlmpc5.calc_initial_guess(ref,traj,5,state,dt); 
    toc

    % AUV6 DLMPC
    tic
    [u6, X6_pred] = dlmpc6.calc_control(ref,6,traj,state,state_pred,u6_0,dt);
    state_pred(:,:,6) = X6_pred;
    u6_actual = u6(1:nu,1);
    U6(:,i) = u6_actual;
    auv6.advance(u6_actual, disturbance, dt);
    X6(:,i+1) = auv6.X;
    state(:,6) = auv6.X;
    u6_0 = dlmpc6.calc_initial_guess(ref,traj,6,state,dt); 
    toc
    
    t = t + dt;
end

tEnd = toc(tStart);

%%
figure(4);

plot(X1(1,end), X1(2,end),'*','color',blue,'LineWidth',2);
hold on;
plot(X2(1,end), X2(2,end),'*','color',lgreen,'LineWidth',2);
hold on;
plot(X3(1,end), X3(2,end),'*','color',red,'LineWidth',2);
hold on;
plot(X4(1,end), X4(2,end),'*','color',yellow,'LineWidth',2);
hold on;
plot(X5(1,end), X5(2,end),'*','color',green,'LineWidth',2);
hold on;
plot(X6(1,end), X6(2,end),'*','color',purple,'LineWidth',2);
hold on;

b3 = plot([X1(1,1),X2(1,1)],[X1(2,1),X2(2,1)],'--','color',cblue);
hold on;
plot([X1(1,1),X3(1,1)],[X1(2,1),X3(2,1)],'--','color',cblue);
hold on;
plot([X1(1,1),X4(1,1)],[X1(2,1),X4(2,1)],'--','color',cblue);
hold on;
plot([X1(1,1),X5(1,1)],[X1(2,1),X5(2,1)],'--','color',cblue);
hold on;
plot([X1(1,1),X6(1,1)],[X1(2,1),X6(2,1)],'--','color',cblue);
hold on;
plot([X2(1,1),X3(1,1)],[X2(2,1),X3(2,1)],'--','color',cblue);
hold on;
plot([X3(1,1),X4(1,1)],[X3(2,1),X4(2,1)],'--','color',cblue);
hold on;
plot([X4(1,1),X5(1,1)],[X4(2,1),X5(2,1)],'--','color',cblue);
hold on;
plot([X5(1,1),X6(1,1)],[X5(2,1),X6(2,1)],'--','color',cblue);
hold on;
plot([X2(1,1),X6(1,1)],[X2(2,1),X6(2,1)],'--','color',cblue);
hold on;

b4 = plot([X1(1,end),X2(1,end)],[X1(2,end),X2(2,end)],'--','color',pink);
hold on;
plot([X1(1,end),X3(1,end)],[X1(2,end),X3(2,end)],'--','color',pink);
hold on;
plot([X1(1,end),X4(1,end)],[X1(2,end),X4(2,end)],'--','color',pink);
hold on;
plot([X1(1,end),X5(1,end)],[X1(2,end),X5(2,end)],'--','color',pink);
hold on;
plot([X1(1,end),X6(1,end)],[X1(2,end),X6(2,end)],'--','color',pink);
hold on;
plot([X2(1,end),X3(1,end)],[X2(2,end),X3(2,end)],'--','color',pink);
hold on;
plot([X3(1,end),X4(1,end)],[X3(2,end),X4(2,end)],'--','color',pink);
hold on;
plot([X4(1,end),X5(1,end)],[X4(2,end),X5(2,end)],'--','color',pink);
hold on;
plot([X5(1,end),X6(1,end)],[X5(2,end),X6(2,end)],'--','color',pink);
hold on;
plot([X2(1,end),X6(1,end)],[X2(2,end),X6(2,end)],'--','color',pink);
hold on;

m1=plot(X1(1,:), X1(2,:),'color',blue,'LineWidth',2);
hold on;
m2=plot(X2(1,:), X2(2,:),'color',lgreen,'LineWidth',2);
hold on;
m3=plot(X3(1,:), X3(2,:),'color',red,'LineWidth',2);
hold on;
m4=plot(X4(1,:), X4(2,:),'color',yellow,'LineWidth',2);
hold on;
m5=plot(X5(1,:), X5(2,:),'color',green,'LineWidth',2);
hold on;
m6=plot(X6(1,:), X6(2,:),'color',purple,'LineWidth',2);
hold on;

legend([m1,m2,m3,m4,m5,m6,b3,b4],...
{'RGMPC UUV1','RGMPC UUV2', 'RGMPC UUV3', 'RGMPC UUV4','RGMPC UUV5','RGMPC UUV6','Initial formation','Final formation'},'NumColumns',4,'location','northeast');
legend('boxoff');
axis square
set(gca,'FontSize',12);
xlabel('\boldmath{$x[m]$}', 'Interpreter', 'latex','FontSize',14);
ylabel('\boldmath{$y[m]$}', 'Interpreter', 'latex','FontSize',14);
set(gcf,'unit','normalized','position', [0.2,0.2,0.6,0.6]);

%%
figure(5);
t = tiledlayout(3,1);
nexttile
m1=plot(0.1:0.1:Tstep/10,U1(1,:),'color',blue,'LineWidth',2);
hold on;
m2=plot(0.1:0.1:Tstep/10,U2(1,:),'color',lgreen,'LineWidth',2);
hold on;
m3=plot(0.1:0.1:Tstep/10,U3(1,:),'color',red,'LineWidth',2);
hold on;
m4=plot(0.1:0.1:Tstep/10,U4(1,:),'color',yellow,'LineWidth',2);
hold on;
m5=plot(0.1:0.1:Tstep/10,U5(1,:),'color',green,'LineWidth',2);
hold on;
m6=plot(0.1:0.1:Tstep/10,U6(1,:),'color',purple,'LineWidth',2);
xlabel('\bf{time}\boldmath{$[s]$}', 'Interpreter', 'latex','FontSize',14);
ylabel('\textbf{control input} \boldmath{$\tau^{x}[N]$}', 'Interpreter', 'latex','FontSize',14);
ylim([-1000 1000]);
legend([m1,m2,m3,m4,m5,m6],...
{'RGMPC UUV1','RGMPC UUV2', 'RGMPC UUV3', 'RGMPC UUV4','RGMPC UUV5','RGMPC UUV6'},'NumColumns',3);
legend('boxoff');
set(gca,'FontSize',12);
nexttile
m1=plot(0.1:0.1:Tstep/10,U1(2,:),'color',blue,'LineWidth',2);
hold on;
m2=plot(0.1:0.1:Tstep/10,U2(2,:),'color',lgreen,'LineWidth',2);
hold on;
m3=plot(0.1:0.1:Tstep/10,U3(2,:),'color',red,'LineWidth',2);
hold on;
m4=plot(0.1:0.1:Tstep/10,U4(2,:),'color',yellow,'LineWidth',2);
hold on;
m5=plot(0.1:0.1:Tstep/10,U5(2,:),'color',green,'LineWidth',2);
hold on;
m6=plot(0.1:0.1:Tstep/10,U6(2,:),'color',purple,'LineWidth',2);
xlabel('\bf{time}\boldmath{$[s]$}', 'Interpreter', 'latex','FontSize',14);
ylabel('\textbf{control input} \boldmath{$\tau^{y}[N]$}', 'Interpreter', 'latex','FontSize',14);
ylim([-1000 1000]);
set(gca,'FontSize',12);
nexttile
m1=plot(0.1:0.1:Tstep/10,U1(3,:),'color',blue,'LineWidth',2);
hold on;
m2=plot(0.1:0.1:Tstep/10,U2(3,:),'color',lgreen,'LineWidth',2);
hold on;
m3=plot(0.1:0.1:Tstep/10,U3(3,:),'color',red,'LineWidth',2);
hold on;
m4=plot(0.1:0.1:Tstep/10,U4(3,:),'color',yellow,'LineWidth',2);
hold on;
m5=plot(0.1:0.1:Tstep/10,U5(3,:),'color',green,'LineWidth',2);
hold on;
m6=plot(0.1:0.1:Tstep/10,U6(3,:),'color',purple,'LineWidth',2);
xlabel('\bf{time}\boldmath{$[s]$}', 'Interpreter', 'latex','FontSize',14);
ylabel('\textbf{control input} \boldmath{$\tau^{r}[N{\cdot}m]$}', 'Interpreter', 'latex','FontSize',14);
set(gca,'FontSize',12);
ylim([-1000 1000]);
set(gcf,'unit','normalized','position', [0.2,0.2,0.5,0.5]);
