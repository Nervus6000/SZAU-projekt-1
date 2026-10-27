%% SZAU - Zadanie 9 - Punkt 1
% Model nieliniowy, linearyzacja, porownanie modeli, analityczny MPC,
% test zbieznosci Eulera. Bez toolboxow sterowania i optymalizacji.
% Wymuszenia od -75% do +100% (bez +125% i +150%).
clear; clc; close all;

%% KONFIGURACJA
outDir = fullfile(pwd,'Wyniki_Zadanie9_Punkt1');
if ~exist(outDir,'dir'), mkdir(outDir); end
p.C1=0.45; p.C2=0.75; p.a1=19; p.a2=18;
p.u0=75; p.d0=10; p.tau=200;
p.h10=((p.u0+p.d0)/p.a1)^2;
p.h20=((p.u0+p.d0)/p.a2)^2;
x0=[p.h10;p.h20];

% Jakobian modelu ciaglego w punkcie pracy
A=[-p.a1/(6*p.C1*p.h10^(5/2)), 0; ...
    p.a1/(4*p.C2*p.h20*sqrt(p.h10)), ...
   -p.a2/(4*p.C2*p.h20^(3/2))];
B=[1/(3*p.C1*p.h10^2);0]; E=B;
fprintf('\n=== PUNKT PRACY I LINEARYZACJA ===\n');
fprintf('h10 = %.6f cm, h20 = %.6f cm\n',p.h10,p.h20);
disp('A =');disp(A);disp('B = E =');disp(B);
disp('Wartosci wlasne A =');disp(eig(A));

%% A. POROWNANIE MODELU NIELINIOWEGO I LINIOWEGO
Tp=1; Tend=12000; t=0:Tp:Tend; K=numel(t);
stepTime=100; delaySteps=round(p.tau/Tp);
changes=[-75 -50 -25 -10 10 25 50 75 100];
nChanges=numel(changes); colors=lines(nChanges);
% Wyniki: rodzaj, procent, RMSE h1,h2, MAX h1,h2, h1,h2 koniec
results=zeros(2*nChanges,8); row=0;
for kind=1:2
    Hnl1=zeros(K,nChanges); Hnl2=Hnl1;
    Hlin1=Hnl1; Hlin2=Hnl1;
    for j=1:nChanges
        pct=changes(j); u=p.u0*ones(1,K); d=p.d0*ones(1,K);
        if kind==1
            u(t>=stepTime)=p.u0*(1+pct/100);
        else
            d(t>=stepTime)=p.d0*(1+pct/100);
        end
        xn=zeros(2,K); xl=zeros(2,K); xn(:,1)=x0; xl(:,1)=x0;
        for k=1:K-1
            idx=k-delaySteps;
            if idx>=1, uDelayed=u(idx); else, uDelayed=p.u0; end
            xn(:,k+1)=xn(:,k)+Tp*tankRHS(xn(:,k),uDelayed,d(k),p);
            if any(xn(:,k+1)<=0) || any(~isfinite(xn(:,k+1)))
                error('Niepoprawny poziom: rodzaj=%d, zmiana=%g%%, t=%g s',kind,pct,t(k+1));
            end
            dx=xl(:,k)-x0;
            xl(:,k+1)=xl(:,k)+Tp*(A*dx+B*(uDelayed-p.u0)+E*(d(k)-p.d0));
        end
        Hnl1(:,j)=xn(1,:)'; Hnl2(:,j)=xn(2,:)';
        Hlin1(:,j)=xl(1,:)'; Hlin2(:,j)=xl(2,:)';
        err=xn-xl; row=row+1;
        results(row,:)=[kind,pct,sqrt(mean(err(1,:).^2)),sqrt(mean(err(2,:).^2)), ...
            max(abs(err(1,:))),max(abs(err(2,:))),xn(1,end),xn(2,end)];
    end
    if kind==1, label='sterowanie'; else, label='zaklocenie'; end
    plotModelFamily(t,Hnl1,Hlin1,changes,colors,p.h10, ...
        ['Zmiany wejscia: ' label ' - h1'],'h_1 [cm]', ...
        fullfile(outDir,['Euler_zbiorczy_' label '_h1']));
    plotModelFamily(t,Hnl2,Hlin2,changes,colors,p.h20, ...
        ['Zmiany wejscia: ' label ' - h2'],'h_2 [cm]', ...
        fullfile(outDir,['Euler_zbiorczy_' label '_h2']));
end
R=array2table(results,'VariableNames',{'wejscie_1u_2d','zmiana_procent', ...
    'RMSE_h1_cm','RMSE_h2_cm','MAX_h1_cm','MAX_h2_cm', ...
    'h1_NL_koniec_cm','h2_NL_koniec_cm'});
fprintf('\n=== A. POROWNANIE MODELI ===\n');disp(R);
writetable(R,fullfile(outDir,'A_porownanie_modeli.csv'));
fig=figure('Color','w','Position',[100 100 1050 540]);
plot(changes,results(1:nChanges,4),'-ob','LineWidth',1.7,'MarkerSize',7);hold on;
plot(changes,results(nChanges+1:end,4),'--sr','LineWidth',1.7,'MarkerSize',7);
xlabel('Zmiana wejscia [%]');ylabel('RMSE h_2 [cm]');
title('Wplyw amplitudy wymuszenia na blad linearyzacji');
legend('Sterowanie F_{1in}','Zaklocenie F_D','Location','northwest');
grid on;saveFigure(fig,fullfile(outDir,'Euler_RMSE_duze_zmiany'));

%% B. ANALITYCZNY MPC - model predykcyjny ZOH, obiekt Euler 1 s
% Parametry wyjsciowe - uczciwy punkt odniesienia, nie 'optymalne' bez strojenia
Ts=20; Np=100; Nu=15; lambda=0.5;
uMin=0; uMax=150;
assert(mod(p.tau,Ts)==0,'Opóznienie musi byc wielokrotnoscia Ts');
assert(mod(Ts,Tp)==0,'Ts musi byc wielokrotnoscia Tp');
L=round(p.tau/Ts); nSub=round(Ts/Tp);
% Dokladna dyskretyzacja liniowego modelu ciaglego dla stalych wejsc
M=expm([A B E;zeros(2,4)]*Ts);
Ad=M(1:2,1:2); Bd=M(1:2,3); Ed=M(1:2,4);

% Stan rozszerzony z kolejka opoznienia i poprzednim sterowaniem:
% z=[delta_h1;delta_h2;u(k-L)-u0;...;u(k-1)-u0;u(k-1)-u0]
% Sterowanie decyzyjne: Delta u(k)=u(k)-u(k-1).
nz=2+L+1;
Aa=zeros(nz); Ba=zeros(nz,1); Ea=zeros(nz,1); Ca=zeros(1,nz);
Aa(1:2,1:2)=Ad; Aa(1:2,3)=Bd;
for j=1:L-1, Aa(2+j,3+j)=1; end
Aa(2+L,nz)=1;  % nowy koniec kolejki = poprzednie u
Aa(nz,nz)=1;   % poprzednie u pozostaje, jesli Delta u=0
Ba(2+L)=1; Ba(nz)=1;
Ea(1:2)=Ed; Ca(2)=1;

% Macierze predykcji Y = F*z + G*DeltaU + D*delta_d
F=zeros(Np,nz); G=zeros(Np,Nu); D=zeros(Np,1);
Aj=eye(nz); disturbAccum=zeros(nz,1);
for j=1:Np
    Aj=Aj*Aa;
    F(j,:)=Ca*Aj;
    disturbAccum=Aa*disturbAccum+Ea;
    D(j)=Ca*disturbAccum;
    for m=1:min(j,Nu)
        G(j,m)=Ca*(Aa^(j-m))*Ba;
    end
end
H=G'*G+lambda*eye(Nu);
% Cholesky - jedna faktoryzacja dla calej symulacji
Hchol=chol(H,'lower');

Tsim=20000; tM=0:Ts:Tsim; Km=numel(tM);
ref=p.h20*ones(1,Km);
ref(tM>=2500)=p.h20+5;
ref(tM>=9000)=p.h20-4;
ref(tM>=15000)=p.h20+3;
dM=p.d0*ones(1,Km);
dM(tM>=6000)=p.d0+5;
dM(tM>=12500)=p.d0-5;

xM=zeros(2,Km);xM(:,1)=x0;
uCmd=p.u0*ones(1,Km); uApplied=p.u0*ones(1,Km);
queue=zeros(L,1); uPrev=0;
for k=1:Km-1
    z=[xM(:,k)-x0;queue;uPrev];
    free=F*z+D*(dM(k)-p.d0);
    desired=(ref(k)-p.h20)*ones(Np,1);
    rhsOpt=G'*(desired-free);
    dU=Hchol'\(Hchol\rhsOpt);
    uNew=min(uMax,max(uMin,p.u0+uPrev+dU(1)));
    uCmd(k)=uNew;
    uApplied(k)=p.u0+queue(1);
    h=xM(:,k);
    for q=1:nSub
        h=h+Tp*tankRHS(h,uApplied(k),dM(k),p);
        if any(h<=0) || any(~isfinite(h))
            error('MPC: niepoprawny poziom przy t=%g s',tM(k));
        end
    end
    xM(:,k+1)=h;
    queue=[queue(2:end);uNew-p.u0];
    uPrev=uNew-p.u0;
end
uCmd(end)=uCmd(end-1);
uApplied(end)=p.u0+queue(1);
eM=xM(2,:)-ref;
RMSE=sqrt(mean(eM.^2)); IAE=trapz(tM,abs(eM)); ISE=trapz(tM,eM.^2);
TVu=sum(abs(diff(uCmd)));
fprintf('\n=== B. ANALITYCZNY MPC (ZOH, Ts=%g, Np=%d, Nu=%d, lambda=%g) ===\n',Ts,Np,Nu,lambda);
fprintf('RMSE = %.6f cm\nIAE = %.6f cm*s\nISE = %.6f cm^2*s\nTVu = %.6f cm^3/s\n',RMSE,IAE,ISE,TVu);
TM=table(tM(:),xM(1,:)',xM(2,:)',ref(:),uCmd(:),uApplied(:),dM(:), ...
    'VariableNames',{'czas_s','h1_cm','h2_cm','h2_zadane_cm','u_zadane','u_opoznione','zaklocenie'});
writetable(TM,fullfile(outDir,'B_MPC_przebiegi.csv'));
metrics=table(Ts,Np,Nu,lambda,RMSE,IAE,ISE,TVu);
writetable(metrics,fullfile(outDir,'B_MPC_wskazniki.csv'));
fig=figure('Color','w','Position',[80 80 1120 550]);
plot(tM,xM(1,:),'b','LineWidth',1.5);hold on;
plot(tM,xM(2,:),'r','LineWidth',1.5);
xlabel('Czas [s]');ylabel('Poziom [cm]');title('MPC - poziomy zbiornikow');
legend('h_1','h_2','Location','best');grid on;
saveFigure(fig,fullfile(outDir,'Euler_MPC_poziomy'));
fig=figure('Color','w','Position',[80 60 1120 850]);
subplot(3,1,1);
plot(tM,xM(2,:),'b','LineWidth',1.7);hold on;
stairs(tM,ref,'k--','LineWidth',1.2);
ylabel('h_2 [cm]');title('Analityczny MPC - obiekt nieliniowy (Euler), predykcja ZOH');
legend('h_2','Wartosc zadana','Location','best');grid on;
subplot(3,1,2);
stairs(tM,uCmd,'m','LineWidth',1.4);hold on;
stairs(tM,uApplied,'g--','LineWidth',1.3);
ylabel('Przeplyw [cm^3/s]');legend('F_{1in} zadane','F_1 opoznione','Location','best');grid on;
subplot(3,1,3);
stairs(tM,dM,'Color',[0.85 0.33 0.1],'LineWidth',1.5);
ylabel('F_D [cm^3/s]');xlabel('Czas [s]');grid on;
saveFigure(fig,fullfile(outDir,'Euler_MPC_regulacja'));

%% C. TEST ZBIEZNOSCI EULERA: -75%, +50%, +100%
testChanges=[-75 50 100];
testSteps=[1 0.5 0.25];
convRes=zeros(numel(testChanges),11);
for j=1:numel(testChanges)
    pct=testChanges(j);
    curves=cell(1,3);
    for q=1:3
        dt=testSteps(q); tc=0:dt:Tend; Nc=numel(tc);
        hc=zeros(2,Nc);hc(:,1)=x0;
        uStep=p.u0*(1+pct/100);
        for k=1:Nc-1
            % Uwaga: skok zadania sterowania w t=100 s, fizyczny w t=300 s
            uNow=p.u0;
            if tc(k)>=stepTime+p.tau, uNow=uStep; end
            hc(:,k+1)=hc(:,k)+dt*tankRHS(hc(:,k),uNow,p.d0,p);
            if any(hc(:,k+1)<=0) || any(~isfinite(hc(:,k+1)))
                error('Test zbieznosci: niepoprawny poziom dla %g%% i dt=%g',pct,dt);
            end
        end
        curves{q}=hc;
    end
    h1=curves{1};h05=curves{2}(:,1:2:end);h025=curves{3}(:,1:4:end);
    e1=h1-h05; e2=h05-h025;
    r1=sqrt(mean(e1.^2,2));r2=sqrt(mean(e2.^2,2));
    m1=max(abs(e1),[],2);m2=max(abs(e2),[],2);
    convRes(j,:)=[pct,r1(1),r2(1),r1(2),r2(2), ...
        m1(1),m2(1),m1(2),m2(2),r1(1)/r2(1),r1(2)/r2(2)];
    fprintf('\nTest zbieznosci %+g%%: RMSE h2=%.8f / %.8f, MAX h2=%.8f, iloraz=%.4f\n', ...
        pct,r1(2),r2(2),m1(2),r1(2)/r2(2));
    fig=figure('Color','w','Position',[80 60 1080 840]);
    for level=1:2
        subplot(3,1,level);
        plot(t,h1(level,:),'b','LineWidth',1.3);hold on;
        plot(t,h05(level,:),'r--','LineWidth',1.1);
        plot(t,h025(level,:),'k:','LineWidth',1.1);
        ylabel(sprintf('h_%d [cm]',level));grid on;
        legend('Euler 1 s','Euler 0.5 s','Euler 0.25 s','Location','best');
        if level==1,title(sprintf('Test zbieznosci Eulera - sterowanie %+g%%',pct));end
    end
    subplot(3,1,3);
    plot(t,e1(2,:),'b','LineWidth',1.2);hold on;
    plot(t,e2(2,:),'r--','LineWidth',1.2);
    xlabel('Czas [s]');ylabel('Roznica h_2 [cm]');grid on;
    legend('Euler 1 s - 0.5 s','Euler 0.5 s - 0.25 s','Location','best');
    if pct<0, name=sprintf('minus_%g',abs(pct));else,name=sprintf('plus_%g',pct);end
    saveFigure(fig,fullfile(outDir,['Euler_zbieznosc_' name '_procent']));
end
C=array2table(convRes,'VariableNames',{'zmiana_procent','RMSE_h1_1_05', ...
    'RMSE_h1_05_025','RMSE_h2_1_05','RMSE_h2_05_025', ...
    'MAX_h1_1_05','MAX_h1_05_025','MAX_h2_1_05', ...
    'MAX_h2_05_025','iloraz_RMSE_h1','iloraz_RMSE_h2'});
fprintf('\n=== C. ZBIORCZY TEST ZBIEZNOSCI ===\n');disp(C);
writetable(C,fullfile(outDir,'C_test_zbieznosci_Eulera.csv'));

fprintf('\n=== KONIEC ===\n');
fprintf('Wyniki: %s\n',outDir);
fprintf('Wykresy: 10 PNG i 10 FIG; tabele CSV.\n');
fprintf('Zakres zmian [%%]: ');fprintf('%g ',changes);fprintf('\n');

%% FUNKCJE LOKALNE
function dh=tankRHS(h,u,d,p)
    if any(h<=0), error('Poziomy h1 i h2 musza byc dodatnie'); end
    dh=[(u+d-p.a1*sqrt(h(1)))/(3*p.C1*h(1)^2); ...
        (p.a1*sqrt(h(1))-p.a2*sqrt(h(2)))/(2*p.C2*h(2))];
end
function plotModelFamily(t,hn,hl,changes,colors,nominal,ttl,ylab,file)
    fig=figure('Color','w','Position',[70 80 1350 700]);hold on;
    n=numel(changes);handles=gobjects(2*n,1);labels=cell(2*n,1);
    for j=1:n
        handles(2*j-1)=plot(t,hn(:,j),'Color',colors(j,:),'LineWidth',1.6);
        handles(2*j)=plot(t,hl(:,j),'--','Color',colors(j,:),'LineWidth',1.1);
        labels{2*j-1}=sprintf('NL %+g%%',changes(j));
        labels{2*j}=sprintf('LIN %+g%%',changes(j));
    end
    yline(nominal,':k','Punkt pracy','HandleVisibility','off');
    title(ttl);xlabel('Czas [s]');ylabel(ylab);grid on;
    legend(handles,labels,'Location','eastoutside','NumColumns',2);
    saveFigure(fig,file);
end
function saveFigure(fig,file)
    savefig(fig,[file '.fig']);
    try
        exportgraphics(fig,[file '.png'],'Resolution',170);
    catch
        saveas(fig,[file '.png']);
    end
end
