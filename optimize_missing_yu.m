function [yopt,uopt,sys,stat] = optimize_missing_yu(y,u,ym,um,lambda,opt)
%
% [yopt,uopt,stat] = optimize_missing_yu(y,u,ym,um,lambda,ess,opt)
%
% computes an optimized output yopt and input uopt from the regularized 
% nuclear norm optimization
% 
%     minimize |[U;Y]|_* + lambda/2 * |y-y_meas|_F^2.
%
% Missing inputs and outputs are indicated by um and ym. Only the measured
% outputs are included in the second term of the optimization objective.
% The missing inputs and outputs are treated as free variables.  
%
% INPUT
% y         p x N real matrix for the outputs, p is the number of outputs,
%           N is the number of data points
% u         m x N real matrix for the inputs, m is the number of inputs
% ym        v x 2 real matrix indicating the missing outputs, each row i
%           indicates y(ym(i,1),ym(i,2)) is missing
% um        w x 2 real matrix indicating the missing inputs, each row i 
%           indicates u(um(i,1),um(i,2)) is missing 
% lambda    regularization parameter, if lambda == Inf, the optimization 
%           reduces to minimize |[U;Y]|_*
% opt       optional structure array for specifying parameters of the 
%           alternating direction methods of multiplier (ADMM) algorithm
% opt.order optional parameter to specify a model order
%    
% OUTPUT
% yopt      optimized output y
% uopt      optimized input u
% sys       state-space model with the following fields
%           -- 'A', 'B', 'C', 'D', the state-space matrices
%           -- 'x1', the initial state x(1)
%           -- 'n', the model order
% stat      structure array with statistics of the ADMM algorithm
%
% The regularization parameter lambda can be entered as a vector with all 
% entries < Inf.  If it is a vector with L entries, L problem instances of 
%     minimize |[U;Y]|_* + lambda(ii)/2 * |y-y_meas|_F^2
% are solved.  The results are stored in struct yopt{ii} and uopt{ii}.


% Problem dimensions
[m,N] = size(u);
p = size(y,1);
if size(y,2) ~= N
    error('Input u and output y must have the same number of data points');
end
v = size(ym,1);
w = size(um,1);

% r: block row dimension in Hankel matrices Y and U
r = 15;
if p == 1, r = 30; end

% Dimensions of block Hankel matrices:  U is rm x M, Y is rp x M
M = N+1-r;

% Optimization variables are
%  [y(1);...;y(N);u(um(1));...;u(um(w))] if lambda ~= INF,
%  [y(ym(1));...;y(ym(v));u(um(1));...;u(um(w))] if lambda == INF

Au = hankel_map(r,M,m,1);
Aum = Au(:,(um(:,2)-1)*m+um(:,1));
Ay = hankel_map(r,M,p,1);
Aym = Ay(:,(ym(:,2)-1)*p+ym(:,1));

if lambda == Inf
    A = @(x) [reshape(Aum*x(v+1:v+w),r*m,M); reshape(Aym*x(1:v),r*p,M)];
    Aa = @(UY) [Aym'*reshape(UY(r*m+1:r*m+r*p,:),r*p*M,1); Aum'*reshape(UY(1:r*m,:),r*m*M,1)];
    C = sparse(zeros(v+w));
    MM = [Aym'*Aym, sparse(v,w); sparse(w,v), Aum'*Aum];
else
    A = @(x) [reshape(Aum*x(p*N+1:p*N+w),r*m,M); hankel_blk(x(1:p*N),r,M,p,1)];
    Aa = @(UY) [hankel_blk_adj(UY(r*m+1:r*m+r*p,:),r,M,p,1);Aum'*reshape(UY(1:r*m,:),r*m*M,1)];
    Cdiag = [ones(p*N,1);zeros(w,1)];
    for i = 1:v
        Cdiag((ym(i,2)-1)*p+ym(i,1)) = 0;
    end
    C = sparse(diag(Cdiag));
    MM = [Ay'*Ay, sparse(p*N,w); sparse(w,p*N), Aum'*Aum];
end
B = [-hankel_blk(u,r,M,m,1);-hankel_blk(y,r,M,p,1)];
x0 = zeros(size(C,1),1);

for ii = 1:length(lambda) 
    % Solve the regularized nuclear norm optimization
    opt.w = lambda(ii);
    if lambda == Inf, opt.w = 1; end
    [dyu, stat1] = rnna_admm(A, Aa, MM, B, C, x0, opt);

    % Compute the optimized output yopt and input uopt
    if lambda == Inf
        yopt1 = y + sparse(ym(:,1),ym(:,2),dyu(1:v),p,N);
        uopt1 = u + sparse(um(:,1),um(:,2),dyu(v+1:end),m,N);
    else
        yopt1 = y + reshape(dyu(1:p*N),p,N);
        uopt1 = u + sparse(um(:,1),um(:,2),dyu(p*N+1:end),m,N);
    end
    
    % Compute the extended observability matrix from QR and SVD of [U;Z;Y]'
    UY = [hankel_blk(uopt1,r,M,m,1);hankel_blk(yopt1,r,M,p,1)];
    s = floor(r/2);
    R = qr([UY(s*m+1:r*m,:); UY(1:s*m,:); UY(r*m+1:end,:)]',0);
    [Uh,Sh] = svd(R((r-s)*m+1:r*m+s*p,r*m+s*p+1:end)'*R((r-s)*m+1:r*m+s*p,(r-s)*m+1:r*m+s*p),'econ');
    sv = diag(Sh);
    
    % Determine model order
    if exist('opt','var') && isfield(opt,'order')
        n = opt.order;
    else
        logsv = log(sv(1:15));
        n = min(10,find(logsv>(max(logsv)+min(logsv))/2,1,'last'));
    end
    
    % Estimate state-space model
    [sys1] = estimate_ss(Uh(:,1:n), yopt1, uopt1);

    % Store results
    if length(lambda) == 1
        yopt = yopt1;
        uopt = uopt1;
        stat = stat1;
        sys = sys1;
    else
        yopt{ii} = yopt1;
        uopt{ii} = uopt1;
        stat{ii} = stat1;
        sys{ii} = sys1;
    end
end
