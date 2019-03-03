% System identification
% using subspace method and nuclear norm optimization.
% Examples are taken from the DaISy database, 

clear all;

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
    
%% Baseline solution

% Identify model using n4sid
model = n4sid(idata,'best','nk',zeros(1,m),'Cov','None','foc','sta');
norig = size(model.A,1);

% Compute validation FIT score
[YHvorig,FITvorig,X0vorig] = compare(vdata,model,Inf);
    
fprintf('Baseline solution from MATLAB n4sid\n');
fprintf('Baseline: n = %d,  FIT = %.1f\n\n',norig,mean(FITvorig));

ii = 3; % IVM

if ii == 1, type = 'MOESP';
elseif ii == 2, type = 'N4SID';
elseif ii == 3, type = 'IVM';
elseif ii == 4, type = 'CVA';
elseif ii == 5, type = 'none';
elseif ii == 6, type = 'noinstrument';
end
            
% Optimize output y using regularized nuclear norm optimization
time = cputime;
lambda = logspace(-3,3,20);
[YoptT,sys,stat] = optimize_y(y(1:N,:)',u(1:N,:)',type,lambda);
fprintf('Nuclear norm optimization with %s weights, optimization time = %.1f\n', type, cputime-time);

KK = zeros(length(lambda),2);
for jj = 1:length(lambda)
    % Apply n4sid to optimized inputs and outputs
    model = n4sid(iddata(YoptT{jj}',u(1:N,:),1),norig,'nk',zeros(1,m),'Cov','None','foc','sta');
    nopt = size(model.A,1);
    
    % Compute validation FIT score
    [YHvopt,FITvopt,X0vopt] = compare(vdata,model,Inf);
    KK(jj,:) = [nopt,mean(FITvopt)];
end

[KKtemp,I] = max(KK,[],1);
fprintf('Optimized: n = %d,  FIT = %.1f \n',KK(I(2),1),KK(I(2),2));
  