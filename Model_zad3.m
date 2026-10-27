% Numeryczny regulator MPC-SL

clear; close all; clc;

%%PARAMETRY OBIEKTU

p.C1 = 0.45;
p.C2 = 0.75;

p.a1 = 19;
p.a2 = 18;

p.u0 = 75;
p.d0 = 10;

p.tau = 200;

% Czas
p.Tp = 1;          % krok Eulera
p.Ts = 20;         % okres regulatora
p.Tend = 20000;

% MPC
p.Np = 100;
p.Nu = 8;
p.lambda = 5;

% Ograniczenia
p.umin = 0;
p.umax = 150;
p.dumax = 3;

% Punkt pracy
p.h10 = ((p.u0+p.d0)/p.a1)^2;
p.h20 = ((p.u0+p.d0)/p.a2)^2;

% Liczba probek opoznienia
p.L = round(p.tau/p.Ts);

assert(abs(p.L*p.Ts-p.tau)<1e-10);
assert(abs(round(p.Ts/p.Tp)*p.Tp-p.Ts)<1e-10);
assert(p.L>=1);

fprintf('Punkt pracy:\n');
fprintf('h1 = %.6f cm\n',p.h10);
fprintf('h2 = %.6f cm\n',p.h20);

%% FOLDER WYNIKOW

outDir = ...
'C:\Users\Julian\Desktop\SZAU_Matlab\Zdje';

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% SCENARIUSZ SYMULACJI

T = (0:p.Ts:p.Tend)';
N = length(T);

% Wartosc zadana
ref = p.h20*ones(N,1);

ref(T>=2500) = p.h20+5;
ref(T>=9000) = p.h20-4;
ref(T>=15000) = p.h20+3;

% Zaklocenie
d = p.d0*ones(N,1);

d(T>=6000) = p.d0+5;
d(T>=12500) = p.d0-5;

%% MODEL LINIOWY W PUNKCIE PRACY

x0 = [p.h10;p.h20];

[A0,B0,E0,~] = localModel(x0,p.u0,p.d0,p);

disp('Macierz A:');
disp(A0);

disp('Macierz B:');
disp(B0);

disp('Macierz E:');
disp(E0);

%% SYMULACJA MPC-SL

% Oba regulatory steruja tym samym obiektem nieliniowym

fprintf('\nUruchamianie MPC analitycznego...\n');

R_MPC = runScenario( ...
    p,T,ref,d,'analytic',A0,B0,E0);

fprintf('\nUruchamianie MPC-SL numerycznego...\n');

R_SL = runScenario( ...
    p,T,ref,d,'SL',A0,B0,E0);

%% WYNIKI REGULACJI

fprintf('\nWYNIKI MPC ANALITYCZNEGO\n');

fprintf('RMSE = %.4f cm\n',R_MPC.RMSE);
fprintf('IAE  = %.4f cm*s\n',R_MPC.IAE);
fprintf('ISE  = %.4f cm^2*s\n',R_MPC.ISE);
fprintf('Nieudane optymalizacje = %d\n',R_MPC.failed);

fprintf('\nWYNIKI MPC-SL NUMERYCZNEGO\n');

fprintf('RMSE = %.4f cm\n',R_SL.RMSE);
fprintf('IAE  = %.4f cm*s\n',R_SL.IAE);
fprintf('ISE  = %.4f cm^2*s\n',R_SL.ISE);
fprintf('Nieudane optymalizacje = %d\n',R_SL.failed);

%% WYKRES REGULACJI

f1 = figure('Color','w', ...
    'Position',[100 100 1150 800]);

subplot(3,1,1);

plot(T,R_MPC.h(2,:),'b','LineWidth',1.6);
hold on;

plot(T,R_SL.h(2,:),'r--','LineWidth',1.5);

plot(T,ref,'k:','LineWidth',1.5);

grid on;

ylabel('h_2 [cm]');

title('Porownanie MPC analitycznego i MPC-SL');

legend('MPC analityczny','MPC-SL numeryczny', ...
    'Wartosc zadana','Location','best');

subplot(3,1,2);

plot(T,R_MPC.u,'b','LineWidth',1.4);
hold on;

plot(T,R_SL.u,'r--','LineWidth',1.4);

grid on;

ylabel('F_{1in} [cm^3/s]');

legend('MPC analityczny','MPC-SL numeryczny', ...
    'Location','best');

subplot(3,1,3);

plot(T,d,'Color',[0.85 0.3 0.1], ...
    'LineWidth',1.5);

grid on;

xlabel('Czas [s]');
ylabel('F_D [cm^3/s]');

saveChart(f1,outDir,'P3_porownanie_regulatorow');

f2 = figure('Color','w', ...
    'Position',[100 100 1100 700]);

subplot(2,1,1);

plot(T,R_MPC.h(1,:),'b','LineWidth',1.5);
hold on;

plot(T,R_SL.h(1,:),'r--','LineWidth',1.5);

yline(p.h10,'k:');

grid on;

ylabel('h_1 [cm]');

title('Poziom pierwszego zbiornika');

legend('MPC analityczny','MPC-SL numeryczny', ...
    'Punkt pracy','Location','best');

subplot(2,1,2);

plot(T,R_MPC.h(2,:),'b','LineWidth',1.5);
hold on;

plot(T,R_SL.h(2,:),'r--','LineWidth',1.5);

plot(T,ref,'k:','LineWidth',1.4);

grid on;

xlabel('Czas [s]');
ylabel('h_2 [cm]');

title('Poziom drugiego zbiornika');

legend('MPC analityczny','MPC-SL numeryczny', ...
    'Wartosc zadana','Location','best');

saveChart(f2,outDir,'P3_porownanie_poziomow');

%% WYKRES BLEDOW REGULACJI

f3 = figure('Color','w', ...
    'Position',[100 100 1100 550]);

plot(T,R_MPC.err,'b','LineWidth',1.5);
hold on;

plot(T,R_SL.err,'r--','LineWidth',1.5);

yline(0,'k:');

grid on;

xlabel('Czas [s]');
ylabel('Blad regulacji [cm]');

title('Porownanie bledow regulacji');

legend('MPC analityczny','MPC-SL numeryczny', ...
    'Location','best');

saveChart(f3,outDir,'P3_porownanie_bledow');

wyniki = table( ...
    T, ...
    R_MPC.h(1,:)', ...
    R_MPC.h(2,:)', ...
    R_SL.h(1,:)', ...
    R_SL.h(2,:)', ...
    ref, ...
    R_MPC.u, ...
    R_SL.u, ...
    R_MPC.uApplied, ...
    R_SL.uApplied, ...
    d, ...
    R_MPC.err, ...
    R_SL.err, ...
    'VariableNames', ...
    {'czas_s', ...
     'h1_MPC_cm','h2_MPC_cm', ...
     'h1_SL_cm','h2_SL_cm', ...
     'zadane_cm', ...
     'u_MPC_cm3s','u_SL_cm3s', ...
     'u_MPC_opoznione_cm3s', ...
     'u_SL_opoznione_cm3s', ...
     'FD_cm3s', ...
     'blad_MPC_cm','blad_SL_cm'});

writetable(wyniki,fullfile(outDir, ...
    'P3_wyniki.csv'));

% Zestawienie wskaznikow jakosci

model = {'MPC analityczny';'MPC-SL numeryczny'};

RMSE = [R_MPC.RMSE;R_SL.RMSE];
IAE = [R_MPC.IAE;R_SL.IAE];
ISE = [R_MPC.ISE;R_SL.ISE];

failed = [R_MPC.failed;R_SL.failed];

jakosc = table(model,RMSE,IAE,ISE,failed);

disp(jakosc);

writetable(jakosc,fullfile(outDir, ...
    'P3_jakosc.csv'));

fprintf('\nZapisano wyniki w folderze:\n%s\n',outDir);

function R = runScenario(p,T,ref,d,controllerType,A0,B0,E0)

N = length(T);
L = p.L;

nx = 2+L+1;

h = zeros(2,N);

h(:,1) = [p.h10;p.h20];

% Sterowanie
u = p.u0*ones(N,1);
uApplied = p.u0*ones(N,1);

% Kolejka opoznienia
queue = p.u0*ones(L,1);

uPrev = p.u0;

exitflag = ones(N-1,1);
opt = optimoptions('quadprog', ...
    'Display','off', ...
    'Algorithm','interior-point-convex');

% Macierz sumujaca przyrosty sterowania
S = tril(ones(p.Nu));

Aineq = [S;-S];

lb = -p.dumax*ones(p.Nu,1);
ub =  p.dumax*ones(p.Nu,1);

x0 = [p.h10;p.h20];
fc0 = -A0*x0-B0*p.u0-E0*p.d0;

[Ad0,Bd0,Ed0,cd0] = discretizeModel( ...
    A0,B0,E0,fc0,p.Ts);

[Az0,Bz0,Ez0,cz0] = extendedModel( ...
    Ad0,Bd0,Ed0,cd0,L);

% Macierze predykcji dla stalego modelu
[G0,~,~] = predictionModel( ...
    Az0,Bz0,Ez0,cz0,p.Np,p.Nu);

Kmpc = (G0'*G0+p.lambda*eye(p.Nu))\G0';

%% PETLA REGULACJI

for k = 1:N-1

    hk = h(:,k);

    uk = queue(1);
    dk = d(k);
    if strcmp(controllerType,'SL')

        [Ac,Bc,Ec,fc] = localModel(hk,uk,dk,p);

        %% DYSKRETYZACJA EULERA

        % Dokladna dyskretyzacja ZOH modelu lokalnego
        [Ad,Bd,Ed,cd] = discretizeModel( ...
            Ac,Bc,Ec,fc,p.Ts);

        %% MODEL ROZSZERZONY

        [Az,Bz,Ez,cz] = extendedModel( ...
            Ad,Bd,Ed,cd,L);

    else

        Az = Az0;
        Bz = Bz0;
        Ez = Ez0;
        cz = cz0;

    end

    % Aktualny stan rozszerzony
    z = [hk;queue;uPrev];

    %% MACIERZE PREDYKCJI

    [G,Yfree,~] = predictionModel( ...
        Az,Bz,Ez,cz,p.Np,p.Nu,z,dk);

    %% WARTOSC ZADANA W HORYZONCIE
    target = ref(k)*ones(p.Np,1);

    %% FUNKCJA CELU

    if strcmp(controllerType,'SL')

        H = 2*(G'*G+p.lambda*eye(p.Nu));

        H = (H+H')/2;

        f = 2*G'*(Yfree-target);

        %% OGRANICZENIA STEROWANIA

        bineq = [ ...
            (p.umax-uPrev)*ones(p.Nu,1);
            (uPrev-p.umin)*ones(p.Nu,1)];

        %%  oPTYMALIZACJA NUMERYCZNA

        [DU,~,ef] = quadprog( ...
            H,f,Aineq,bineq, ...
            [],[],lb,ub,[],opt);

        exitflag(k) = ef;

        if ef<=0 || isempty(DU)

            warning( ...
                'Problem quadprog w kroku %d',k);

            DU = zeros(p.Nu,1);

        end

    else

        % Analityczne rozwiazanie MPC bez ograniczen
        DU = Kmpc*(target-Yfree);

    end

    %% NOWE STEROWANIE

    uCmd = uPrev+DU(1);

    % Ograniczenia fizyczne sterowania
    uCmd = min(max(uCmd,p.umin),p.umax);

    % Ograniczenie przyrostu sterowania
    uCmd = min(max(uCmd, ...
        uPrev-p.dumax),uPrev+p.dumax);

    u(k) = uCmd;

    uApplied(k) = queue(1);

    %% OBIEKT NIELINIOWY - EULER 1 s

    steps = round(p.Ts/p.Tp);

    hn = hk;

    for m = 1:steps

        dh = plantRhs( ...
            hn,queue(1),dk,p);

        hn = hn+p.Tp*dh;

        hn = max(hn,1e-6);

    end

    h(:,k+1) = hn;

    %% AKTUALIZACJA OPOZNIENIA

    queue = [queue(2:end);uCmd];

    uPrev = uCmd;

end

%% WYNIKI

u(end) = u(end-1);

uApplied(end) = queue(1);

err = h(2,:)'-ref;

R.h = h;
R.u = u;
R.uApplied = uApplied;
R.err = err;

R.RMSE = sqrt(mean(err.^2));

R.IAE = p.Ts*sum(abs(err(1:end-1)));

R.ISE = p.Ts*sum(err(1:end-1).^2);

R.failed = sum(exitflag<=0);

end

%% MODEL NIELINIOWY

function dh = plantRhs(h,u,d,p)

h = max(h,1e-6);

h1 = h(1);
h2 = h(2);

% Przeplyw z pierwszego do drugiego
F2 = p.a1*sqrt(h1);

% Odplyw z drugiego zbiornika
F3 = p.a2*sqrt(h2);

% Pierwszy zbiornik
dh1 = (u+d-F2)/(3*p.C1*h1^2);

% Drugi zbiornik
dh2 = (F2-F3)/(2*p.C2*h2);

dh = [dh1;dh2];

end

%% SUKCESYWNA LINEARYZACJA

function [A,B,E,fc] = localModel(h,u,d,p)

h = max(h,1e-6);

h1 = h(1);
h2 = h(2);

f0 = plantRhs(h,u,d,p);

%% Macierz A

den1 = 3*p.C1*h1^2;

num1 = u+d-p.a1*sqrt(h1);

A11 = -p.a1/(2*sqrt(h1)*den1) ...
      -2*num1/(3*p.C1*h1^3);

A12 = 0;

A21 = p.a1/(4*p.C2*h2*sqrt(h1));

num2 = p.a1*sqrt(h1)-p.a2*sqrt(h2);

A22 = -p.a2/(4*p.C2*h2^(3/2)) ...
      -num2/(2*p.C2*h2^2);

A = [A11 A12;
     A21 A22];

%% Macierze wejsc

B = [1/den1;0];

E = B;

%% Wyraz afiniczny

fc = f0-A*h-B*u-E*d;

end

%% DYSKRETYZACJA MODELU

function [Ad,Bd,Ed,cd] = discretizeModel(A,B,E,fc,Ts)

% Macierz rozszerzona 5x5
M = [A B E fc;
     zeros(3,5)];

% Dokladna dyskretyzacja ZOH
Md = expm(M*Ts);

Ad = Md(1:2,1:2);

Bd = Md(1:2,3);
Ed = Md(1:2,4);

cd = Md(1:2,5);

end

%% MODEL ROZSZERZONY

function [Az,Bz,Ez,cz] = extendedModel(Ad,Bd,Ed,cd,L)

nx = 2+L+1;

Az = zeros(nx);

Bz = zeros(nx,1);
Ez = zeros(nx,1);
cz = zeros(nx,1);

% Model obiektu
Az(1:2,1:2) = Ad;
Az(1:2,3) = Bd;

Ez(1:2) = Ed;
cz(1:2) = cd;

% Przesuwanie kolejki opoznienia
if L>1
    Az(3:1+L,4:2+L) = eye(L-1);
end

% Nowa komenda
Az(2+L,3+L) = 1;
Az(3+L,3+L) = 1;

Bz(2+L) = 1;
Bz(3+L) = 1;

end

%% MACIERZE PREDYKCJI

function [G,Yfree,F] = predictionModel( ...
    Az,Bz,Ez,cz,Np,Nu,z,dk)

nx = size(Az,1);

G = zeros(Np,Nu);

F = zeros(Np,nx);

sensitiv = zeros(nx,Nu);

Apow = eye(nx);

for j = 1:Np

    Apow = Az*Apow;

    F(j,:) = Apow(2,:);

    sensitiv = Az*sensitiv;

    if j<=Nu
        sensitiv(:,j) = ...
            sensitiv(:,j)+Bz;
    end

    G(j,:) = sensitiv(2,:);

end

if nargin<7
    Yfree = [];
    return;
end

Yfree = zeros(Np,1);

zz = z;

for j = 1:Np

    zz = Az*zz+Ez*dk+cz;

    Yfree(j) = zz(2);

end

end

%% ZAPIS WYKRESOW

function saveChart(f,outDir,name)

savefig(f,fullfile(outDir,[name '.fig']));

exportgraphics(f, ...
    fullfile(outDir,[name '.png']), ...
    'Resolution',200);

end