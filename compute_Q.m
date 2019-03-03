function [Q,d,time] = compute_Q(C,M,method)
%
% [Q,d,time] = compute_Q(C,M,method);
%
% Compute the dense matrix Q that simultaneously diagonalizes C and M
%     Q' C Q = diag(d),    Q' M Q = I - diag(d),
% where C is a symmetric matrix and C+M is a positive definite matrix.
%
% Optional input argument 'method' to specify which method
%     1 for generalized eigenvalue decomposition 
%     2 for Cholesky factorization plus eigenvalue decomposition


% Initialize method
if ~exist('method','var'), method = 1; end

time = cputime;

if method == 1
    %if min(eig(C+M)) < 1E-6,error('C+M must be positive definite'); end
	[Q,D] = eig(full(C),C+M,'chol');
elseif method == 2
	R = chol(C+M); 
    [Q,D] = eig((R')\(C/R)); 
    Q = R\Q; 
end
d = diag(D);

time = cputime - time;

%% Testing script
% n=1000;
% C =diag(rand(n,1));
% M = randn(n);
% M = M*M';
% 
% [Q1,d1,time1] = compute_Q(C,M,1);
% norm(Q1'*C*Q1 - diag(d1),'fro')
% norm(Q1'*M*Q1 - (eye(n)-diag(d1)),'fro')
% 
% [Q2,d2,time2] = compute_Q(C,M,2);
% norm(Q2'*C*Q2 - diag(d2),'fro')
% norm(Q2'*M*Q2 - (eye(n)-diag(d2)),'fro')
