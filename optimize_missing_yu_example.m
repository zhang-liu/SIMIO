% System identification of systems with missing inputs and outputs
% using subspace method and nuclear norm optimization.
% Data is taken from the DaISy database's CD player arm example.

clear all;
rng(1);    

% CD_player_arm example
load CD_player_arm.dat;
u = CD_player_arm(:,1:2);
y = CD_player_arm(:,3:4);
N = 200;
Nt = 600;
    
% Data dimensions
m = size(u,2);
p = size(y,2);
u = u(1:Nt,:);
y = y(1:Nt,:);

fprintf('--- CD Player Arm Example ---\n');
fprintf('Two inputs and two outputs\n');
fprintf('200 data points for identification and 400 data points for validation\n\n');
    
% Inputs and outputs for identification
idata = iddata(y(1:N,:),u(1:N,:),1);

% Inputs and outputs for validation
vdata = iddata(y(N+1:Nt,:),u(N+1:Nt,:),1);
    
%% Apply Matlab N4SID on complete input-output data

% Identify model using n4sid
model = n4sid(idata,'best','nk',zeros(1,m),'Cov','None','foc','sta');
norig = size(model.A,1);

% Compute validation FIT score
[YHv,FITv,X0v] = compare(vdata,model,Inf);

fprintf('Baseline solution using MATLAB n4sid on complete input-output data\n');
fprintf('Baseline: n=%d,  FIT = %.1f\n',norig,mean(FITv));
    
%% Missing inputs and outputs
for missing_percent = 0.2:0.3:0.8
    % Creating missing inputs and outputs at random places
    um_percent = missing_percent;
    ym_percent = missing_percent;
    temp = randperm(m*N)';
    temp = temp(1:ceil(m*N*um_percent));
    um = [mod(temp-1,m)+1, ceil(temp/m)];
    temp = randperm(p*N)';
    temp = temp(1:ceil(p*N*ym_percent));
    ym = [mod(temp-1,p)+1, ceil(temp/p)];
    fprintf('\nPercentage of missing inputs and outputs = %d%% \n',missing_percent*100);
    
    % Set missing inputs and outputs to 0
    Yi = y(1:N,:)';
    for i = 1:size(ym,1)
        Yi(ym(i,1),ym(i,2)) = 0;
    end
    Ui = u(1:N,:)';
    for i = 1:size(um,1)
        Ui(um(i,1),um(i,2)) = 0;
    end       
   
    % Use linear interpolation to fill the missing inputs and outputs
    Yint = zeros(p,N);
    Uint = zeros(m,N);
    t = 1:N;
    ya = ones(p,N);
    for i = 1:size(ym,1)
        ya(ym(i,1),ym(i,2)) = 0;
    end
    ua = ones(m,N);
    for i = 1:size(um,1)
        ua(um(i,1),um(i,2)) = 0;
    end   
    for pp = 1:p
        Yint(pp,:) = interp1(t(ya(pp,:)==1),Yi(pp,ya(pp,:)==1),t,'linear','extrap');
    end
    for mm = 1:m
        Uint(mm,:) = interp1(t(ua(mm,:)==1),Ui(mm,ua(mm,:)==1),t,'linear','extrap');
    end
    
    % Apply n4sid to the linear interpolated data
    model = n4sid(iddata(Yint',Uint',1),norig,'nk',zeros(1,m),'Cov','None','foc','sta');
    n = size(model.A,1);
    [YHv,FITv,X0v] = compare(vdata,model,Inf);
    fprintf('Interpolated: n=%d,  FIT = %.1f\n',n,mean(FITv));
    
    % Optimize output y and input u using regularized nuclear norm optimization
    time = cputime;
    lambda = Inf;
    [Yopt,Uopt,sys,stat] = optimize_missing_yu(Yi,Ui,ym,um,lambda);
    opt_time = cputime-time;
        
    % Apply n4sid to the nuclear norm optimized data
    model = n4sid(iddata(Yopt',Uopt',1),norig,'nk',zeros(1,m),'Cov','None','foc','sta');
    n = size(model.A,1);
    [YHv,FITv,X0v] = compare(vdata,model,Inf);    
    fprintf('Optimized: n = %d,  FIT = %.1f  sol_time = %.1f\n',n,mean(FITv),opt_time);
end 

figure(1);
iiu = 1;
valid = find(Ui(iiu,:) ~= 0);
plot(1:N,u(1:N,iiu),valid,Ui(iiu,valid),'o',1:N,Uopt(iiu,:));
legend('Original','Available','Optimized','location','best');
xlabel('Time t','FontSize',12);
ylabel('Input 1','FontSize',12);
set(gca,'FontSize',12);

figure(2);
iiu = 2;
valid = find(Ui(iiu,:) ~= 0);
plot(1:N,u(1:N,iiu),valid,Ui(iiu,valid),'o',1:N,Uopt(iiu,:));
legend('Original','Available','Optimized','location','best');
xlabel('Time t','FontSize',12);
ylabel('Input 2','FontSize',12);
set(gca,'FontSize',12);

figure(3);
iiy = 1;
valid = find(Yi(iiy,:) ~= 0);
plot(1:N,y(1:N,iiy),valid,Yi(iiy,valid),'o',1:N,Yopt(iiy,:));
legend('Original','Available','Optimized','location','best');
xlabel('Time t','FontSize',12);
ylabel('Output 1','FontSize',12);
set(gca,'FontSize',12);

figure(4);
iiy = 2;
valid = find(Yi(iiy,:) ~= 0);
plot(1:N,y(1:N,iiy),valid,Yi(iiy,valid),'o',1:N,Yopt(iiy,:));
legend('Original','Available','Optimized','location','best');
xlabel('Time t','FontSize',12);
ylabel('Output 2','FontSize',12);
set(gca,'FontSize',12);

