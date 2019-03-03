function [x, stat] = rnna_admm(A, Aa, M, B, C, x0, opt)
% FUNCTION: [x, stat] = rnna_admm(A, Aa, M, B, C, x0, opt);
%
% Implementation of alternating direction method of multipliers (ADMM) for 
% solving the regularized nuclear norm approximation problem.
%
%     minimize |A(x) - B|_* + 1/2 (x - x0)' C (x - x0)
%
% INPUT:
% A           function handle that maps a n x 1 real vector to a p x q
%             real matrix
% Aa          function handle that maps a p x q real matrix to a n x 1 
%             real vector, Aa is the adjoint of A
% M           n x n real positive definite matrix, define by M*x = Aa(A(x))
% B           p x q real matrix
% C           n x n real positive semidefinite matrix
% x0          n x 1 real vector
% opt.maxiter maximum number of iterations, default = 200
% opt.epsabs  absolute solution accuracy tolerance, default = 1E-6
% opt.epsrel  relative solution accuracy tolerance, default = 1E-3
% opt.t       penalty parameter, t > 0
%             default = -1, update t based on primal and dual residuals
% opt.t_mu    coefficient for condition of updating t, default = 10
% opt.t_tau   coefficient for updating t, default = 2  
% opt.w       scalar weight, if specified, the optimization becomes
%                 minimize |A(x) - B|_* + 1/2 w (x - x0)' C (x - x0)
%             default w = 1
% opt.Q       simultaneously diagonalize C and M
% opt.d           Q'*C*Q = diag(d), Q'*M*Q = I - diag(d)
%
% OUTPUT:
% x           solution
% stat.iters  number of iterations
% stat.Ttime  total time in iterations (seconds)
% stat.time   time in seconds at each iteration since start of 1st iter
% stat.rp     primal residual at each iteration
% stat.rd     dual residual at each iteration
% stat.epsp   primal residual tolerance at each iteration
% stat.epsd   dual residual tolerance at each iteration
% stat.t      penalty parameter at each iteration


% Initialize options
if exist('opt','var') && isfield(opt,'maxiter'), maxiter = opt.maxiter; else maxiter = 200; end
if exist('opt','var') && isfield(opt,'epsabs'),  epsabs = opt.epsabs;   else epsabs = 1E-6; end
if exist('opt','var') && isfield(opt,'epsrel'),  epsrel = opt.epsrel;   else epsrel = 1E-3; end
if exist('opt','var') && isfield(opt,'t'),       t = opt.t;             else t = -1;        end
if exist('opt','var') && isfield(opt,'t_mu'),    t_mu = opt.t_mu;       else t_mu = 10;     end
if exist('opt','var') && isfield(opt,'t_tau'),   t_tau = opt.t_tau;     else t_tau = 2;     end
if exist('opt','var') && isfield(opt,'w'),       w = opt.w;             else w = 1;         end
if exist('opt','var') && isfield(opt,'Q') && isfield(opt,'d'), factor = 1; else factor = 0; end

% Determine problem size
[p,q] = size(B);
n = size(x0,1);

% Preprocessing
C = w * C;
Cx0 = C * x0;
normB = norm(B,'fro');

% Initialize statistics variables
stat.rp   = zeros(maxiter,1);
stat.rd   = zeros(maxiter,1);
stat.epsp = zeros(maxiter,1);
stat.epsd = zeros(maxiter,1);
stat.t    = zeros(maxiter,1);
stat.time = zeros(maxiter,1);

%%% Start the ADMM algorithm
stat.Ttime = cputime;

% Initialize variables
x = zeros(n,1);
Y = -B;
Z = zeros(p,q);
if t > 0, t_upd = 0; else t = 1; t_upd = 1; end
t_chg = 1;

% ADMM iterations
for iter = 1:maxiter
    % Update x
    if factor == 1
        x = opt.Q * ( (opt.Q'*(Aa(t*(Y+B)-Z)+Cx0)) ./ ((w-t)*opt.d+t) );
    else
        if t_chg == 1, CtH_R = chol(C + t*M); t_chg = 0; end
        x = CtH_R \ (CtH_R' \ (Aa(t*(Y+B)-Z)+Cx0));
    end
    Ax = A(x);
    
    % Update Y
    Yp = Y;
    [U,S,V] = svd(Ax-B+Z/t);
    Y = U*max(0,S-1/t)*V';
    
    % Update residuals
    rp = Ax - Y - B;
    rd = t*Aa(Yp-Y);
    normrp = norm(rp,'fro');
    normrd = norm(rd);
    
    % Update Z
    Z = Z + t*rp;
    
    % Update tolerance
    epsp = epsrel*max([norm(Ax,'fro'),norm(Y,'fro'),normB]) + sqrt(p*q)*epsabs;
    epsd = (epsrel*norm(Aa(Z)) + sqrt(n)*epsabs)*10;
    
    % Store statistics
    stat.rp(iter) = normrp;
    stat.rd(iter) = normrd;
    stat.epsp(iter) = epsp;
    stat.epsd(iter) = epsd;
    stat.t(iter) = t;
    stat.time(iter) = cputime - stat.Ttime;
    
    % Check stopping criteria
    if normrp <= epsp && normrd <= epsd
        break;
    end
    
    % Update penalty parameter t
    if t_upd == 1
        if (normrp) > t_mu * (normrd)
            t = t_tau * t;
            t_chg = 1;
        elseif (normrd) > t_mu * (normrp)
            t = t / t_tau;
            t_chg = 1;
        end
    end
end

% Store statistics
stat.iters = iter;
stat.Ttime = cputime - stat.Ttime;

stat.rp = stat.rp(1:iter);
stat.rd = stat.rd(1:iter) ;
stat.epsp = stat.epsp(1:iter);
stat.epsd = stat.epsd(1:iter);
stat.t = stat.t(1:iter);
stat.time = stat.time(1:iter);
