function [M,time] = compute_M(r,N,u,L,R,method)
%
% [M,time] = compute_M(r,N,u,L,R,method);
%
% Compute the scaled Gram matrix M, defined as M x = A_adj(A(x)), where 
% A is a linear mapping A(x) = L*Hankel(x)*R and 
%
%             [ H_1   H_2   ...   H_N     ]
%             [ H_2   H_3   ...   H_N+1   ]
% Hankel(x) = [  :     :    ...    :      ]
%             [ H_r   H_r+1 ...   H_N+r-1 ]
%
% is the block Hankel mapping with dimensions of H_i u x 1.
%
% Optional input argument 'method' to specify which method
%     1 for the standard method of forming A and then M = A'*A 
%     2 for the discrete fourier transform (DFT) matrix approach
%     3 for the fast fourier transform approach (default)


% Problem dimension
n = (r+N-1)*u;    
p = size(L,1);
q = size(R,2);
    
% Initialize method
if ~exist('method','var'), method = 3; end

start = cputime;
if method == 1
	H = hankel_map(r,N,u,1);
	A = sparse(p*q,n);
	for ii = 1:n
		temp = L*reshape(H(:,ii),r*u,N)*R;
		A(:,ii) = temp(:);
	end
	M = A'*A;
elseif method == 2
	K = 2*r+2*N-3;
	T = fft(eye(r+N-1),K);
	Tt = fliplr(T);
	E = kron(Tt(:,1:r),eye(u));
	F = kron(Tt,eye(u));
	G = kron(T(:,1:N),ones(u,1));
	M = real(1/K^2 * F' * ((E*L'*L*E').*conj(G*R*R'*G')) * F);
    M = 0.5*(M+M');
elseif method == 3
	K = 2*r+2*N-3;
	% G * R
	GR = kron(fft(R,K),ones(u,1));
	% E * L'
	ELT = zeros(K*u,size(L,1));
	for ii = 1:u
		ELT(ii:u:end,:) = fft([zeros(N-1,size(L,1));flipud(L(:,ii:u:end)')],K);
	end
	% ELGR = (ELT*ELT').*conj(GR*GR')
	ELGR = (ELT*ELT').*conj(GR*GR');
	% MM = 1/K * F' * ELGR
	MM = zeros(n,K*u);
	for ii = 1:u
		temp = ifft(ELGR(ii:u:end,:),K);
		MM(ii:u:end,:) = flipud(temp(1:r+N-1,:));
	end
	% M = 1/K * MM * F
	M = zeros(n,n);
	for ii = 1:u
		temp = ifft(MM(:,ii:u:end)',K);
		M(ii:u:end,:) = flipud(temp(1:r+N-1,:));
	end
	M = real(M');
    M = 0.5*(M+M');
end

time = cputime - start;
    