function H = hankel_map(m, n, p, q)
%
% [H] = hankel_map(m, n, p, q);
%
% Return matrix H that maps a matrix of size (p,q*(m+n-1))
%
%      [ h_1, h_2, ..., h_{m+n-1} ]
%
% stored in column-major order to a block Hankel matrix of size m*p x n*q
%
%      [ h_1   h_2     ...  h_n       ]
%      [ h_2   h_3     ...  h_{n+1}   ]
%      [  :     :            :        ]
%      [ h_m   h_{m+1} ...  h_{m+n-1} ]
%
% stored in column-major order, where h_i is of size p x q.
%   
% INPUT
% m       number of block rows
% n       number of block columns
% p       number of rows of a Hankel block
% q       number of columns of a Hankel block
%
% OUTPUT
% H       matrix of size m*n*p*q x p*q*(m+n-1)

a = m*p;
b = n*q;

Hj = zeros(a,b);
for ii = 0:m-1
    Hj(ii*p+1:(ii+1)*p,:) = reshape((1:p*b)+ii*p*q,p,b);
end
H = sparse(1:a*b,Hj(:),1);
