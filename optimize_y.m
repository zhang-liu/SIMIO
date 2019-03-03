function [yopt,sys,stat] = optimize_y(y,u,type,lambda,opt)
%
% [yopt,sys,stat] = optimize_y(y,u,type,lambda,opt)
%
% computes an optimized output yopt from regularized nuclear norm 
% optimization
% 
%     minimize |W_1*Y*U^perp*Phi^T*W_2|_* + lambda/2 * |y-y_meas|_F^2.
%
% where the weighting matrices W1 and W2 are determined by the input
% variable 'type', which can have values:
%
% 'MOESP', 'N4SID', 'IVM', 'CVA', 'none', 'noinstrument'
%
% When type = 'none', W1 = W2 = I. 
%
% When type = 'noinstrument', W1 = W2 = Phi = I
%
% Note: the first s outputs yopt(:,1:s) used in Yp are not optimized.
%
% INPUT
% y         p x N real matrix for the measured outputs, p is the number
%           of outputs, N is the number of data points
% u         m x N real matrix for the inputs, m is the number of inputs
% opt       optional structure array for specifying parameters of the 
%           alternating direction methods of multiplier (ADMM) algorithm
% opt.order optional parameter to specify a model order
%    
% OUTPUT
% yopt      optimized output y
% sys       state-space model with the following fields
%           -- 'A', 'B', 'C', 'D', the state-space matrices
%           -- 'x1', the initial state x(1)
%           -- 'n', the model order
% stat      structure array with statistics of the ADMM algorithm
%
% The regularization parameter lambda can be entered as a vector.  If it 
% is a vector with L elements, L problem instances of 
%     minimize |W_1*Y*U^perp*Phi^T*W_2|_* + lambda(ii)/2 * |y-y_meas|_F^2
% are solved.  The results are stored in struct yopt{ii},sys{ii},stat{ii}.
% If L > 6, simultaneous diagonalization of C and M is computed.
 
% Dimensions of input and output data
[m,N] = size(u);
p = size(y,1);
if size(y,2) ~= N
    error('Input u and output y must have the same number of data points');
end

% r: block row dimension in Hankel matrices Yf and Uf
% s: block row dimension in Hankel matrices Yp and Up
r = 15;
s = 15;

% Dimensions of block Hankel matrices
% Up: sm x M, Yp: sp x M, Uf: rm x M, Yf: rp x M
M = N+1-r-s;

% Form the Hankel matrices, U = [Up;Uf], Y = [Yp;Yf]
U = hankel_blk(u,s+r,M,m,1);
Y = hankel_blk(y,s+r,M,p,1);
Uf = U(s*m+1:end,:);
Yf = Y(s*p+1:end,:);

% Form the projection matrix
Uperp = eye(M) - Uf'*((Uf*Uf')\Uf);

% Define the instrument variable, Phi = [Yp;Up]
Phi = [Y(1:s*p,:); U(1:s*m,:)];

% Define the weighting matrices W1 and W2
if strcmp(type,'none')
    W1 = eye(r*p);
    W2 = eye(s*(p+m));
elseif strcmp(type,'MOESP')
    W1 = eye(r*p);
    W2 = (Phi*Uperp*Phi')\(Phi*Uperp);
elseif strcmp(type,'N4SID')
    W1 = eye(r*p);
    W2 = (Phi*Uperp*Phi')\Phi;
elseif strcmp(type,'IVM')
    W1 = sqrtm(inv(Yf*Uperp*Yf'));
    W2 = sqrtm(inv(Phi*Phi'));
elseif strcmp(type,'CVA')
    W1 = sqrtm(inv(Yf*Uperp*Yf'));
    W2 = sqrtm(inv(Phi*Uperp*Phi'));
end

% Define the optimization data A, Aa, M, B, C, x0 for rnna_admm.m
% The optimization variable is dy such that yopt = y + [zeros(p,s),dy]

if strcmp(type,'noinstrument')
    UpPW2 = Uperp;
    W2 = speye(size(Uperp,2));
    W1 = eye(r*p);
else
    UpPW2 = Uperp*Phi'*W2;
end

A  = @(x) W1*hankel_blk(x,r,M,p,1)*UpPW2;
Aa = @(Y) hankel_blk_adj(W1'*Y*UpPW2',r,M,p,1);
[MM,M_time] = compute_M(r,M,p,W1,UpPW2);

B = -W1*Yf*UpPW2;
x0 = zeros(p*(N-s),1);
C = speye(p*(N-s));

if length(lambda) > 6
    [opt.Q,opt.d,Q_time] = compute_Q(C,MM);
end

for ii = 1:length(lambda)
    % Solve the regularized nuclear norm optimization
    opt.w = lambda(ii);
    [dy, stat1] = rnna_admm(A, Aa, MM, B, C, x0, opt);
    
    % Compute the optimized outputs yopt
    yopt1 = y + [zeros(p,s),reshape(dy,p,N-s)];
    
    % Optimized singular vlaues
    [Uh,Sh] = svd(A(dy)-B,'econ');
    sv = diag(Sh);
    
    % Determine model order
    if exist('opt','var') && isfield(opt,'order')
        n = opt.order;
    else
        logsv = log(sv(1:15));
        n = min(10,find(logsv>(max(logsv)+min(logsv))/2,1,'last'));
    end
    
    % Estimate state-space model
    [sys1] = estimate_ss(W1\Uh(:,1:n), yopt1(:,s+1:end), u(:,s+1:end));
    sys1.n = n;
    
    % store solutions
    if length(lambda) == 1
        yopt = yopt1;
        stat = stat1;
        sys = sys1;
    else
        yopt{ii} = yopt1;
        stat{ii} = stat1;
        sys{ii} = sys1;
    end
end
