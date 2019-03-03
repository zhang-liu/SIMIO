function [sys] = estimate_ss(Obs, y, u);
%
% Estimate the state space model
%
% INPUT
% Obs      extended observability matrix
% y        output
% u        input 
%
% OUTPUT
% sys.A,B,C,D  state-space model
% sys.x1       initial state x(1)
% sys.n        model order

n = size(Obs,2);
[p,N] = size(y);
m = size(u,1);
             
% estimate A, C
C = Obs(1:p,:);
A = Obs(1:end-p,:)\Obs(p+1:end,:);
Aerr = norm(Obs(1:end-p,:)*A - Obs(p+1:end,:),'fro')/norm(Obs(p+1:end,:),'fro');
               
% stabilize A
[Us,Ts] = schur(A,'complex');
dTs = diag(Ts);
idTs = find(abs(dTs) >= 1.0);
while ~isempty(idTs)
    %disp('Stabilize matrix A');
    dTs(idTs) = dTs(idTs).^-1;
    Ts(1:n+1:end) = dTs;
    A = real(Us*Ts*Us');

    [Us,Ts] = schur(A,'complex');
    dTs = diag(Ts);
    idTs = find(abs(dTs) >= 1.0);
end

% estimate B,D,x1 stored in vector [x1; vec(D); vec(B)]
F1 = zeros(p*N,n);
F1(1:p,:) = C;
for kk = 1:N-1
    F1(kk*p+1:(kk+1)*p,:) = F1((kk-1)*p+1:kk*p,:)*A;
end
F2 = kron(u',eye(p));
F3 = zeros(p*N,n*m);
for kk = 1:N-1
    F3(kk*p+1:end,:) = F3(kk*p+1:end,:) + kron(u(:,1:N-kk)',F1((kk-1)*p+1:kk*p,:));
end
F = [F1,F2,F3];
x = F\y(:);
xerr = norm(F*x-y(:),'fro')/norm(y,'fro');
x1 = x(1:n);
D = reshape(x(n+1:n+p*m),p,m);
B = reshape(x(n+p*m+1:end),n,m);

% store solutions
sys.A = A;
sys.B = B;
sys.C = C;
sys.D = D;
sys.x1 = x1;
sys.n = n;
sys.Aerr = Aerr;
sys.xerr = xerr;
