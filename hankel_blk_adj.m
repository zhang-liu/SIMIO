function x = hankel_blk_adj(H, m, n, p, q)
%
% [x] = hankel_blk_adj(H, m, n, p, q);
%
% Adjoint of the mapping hankel_blk(x,m,n,p,q).  It is the sum of all the 
% block anti-diagonals of matrix H.
%
% INPUT
% H       block Hankel matrix of size (m*p,n*q),
%                [ h_1   h_2     ...  h_n       ]
%            H = [ h_2   h_3     ...  h_{n+1}   ]
%                [  :     :            :        ]
%                [ h_m   h_{m+1} ...  h_{m+n-1} ],
%         where h_i, i=1...(m+n-1), are the Hankel blocks of size (p,q).
% m       number of block rows
% n       number of block columns
% p       number of rows of a Hankel block
% q       number of columns of a Hankel block
%
% OUTPUT
% x       matrix of size (p*q*(m+n-1),1) with the corresponding Hankel
%         blocks stored in column-major order


if size(H,1) ~= m*p || size(H,2) ~= n*q
    error('Error: dimensions of H must equal to (m*p,n*q)');
end

x = zeros(p, q*(m+n-1));
for ii = 1:m
    x(:,(ii-1)*q+1:(ii+n-1)*q) = x(:,(ii-1)*q+1:(ii+n-1)*q) + H((ii-1)*p+1:ii*p,:);
end

x = reshape(x, p*q*(m+n-1), 1);