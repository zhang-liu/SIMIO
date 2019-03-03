function H = hankel_blk(x, m, n, p, q)
%
% [H] = hankel_blk(x, m, n, p, q);
%
% Return the block Hankel matrix H, a matrix of size (m*p,n*q),
%    
%         [ h_1   h_2     ...  h_n       ]
%     H = [ h_2   h_3     ...  h_{n+1}   ]
%         [  :     :            :        ]
%         [ h_m   h_{m+1} ...  h_{m+n-1} ],
%    
% where h_i, i=1...(m+n-1), are the Hankel blocks of size (p,q).
%   
% INPUT
% x       matrix of size (p*q*(m+n-1),1) with the Hankel blocks,
%         [ h_1, h_2, ..., h_{m+n-1} ], stored in column-major order
% m       number of block rows
% n       number of block columns
% p       number of rows of a Hankel block
% q       number of columns of a Hankel block    


if numel(x) ~= p*q*(m+n-1)
    error('Error: number of elements of x must equal to p*q*(m+n-1)');
end

x = reshape(x, p, q*(m+n-1));

H = zeros(m*p,n*q);
for ii = 1:m
    H((ii-1)*p+1:ii*p,:) = x(:,(ii-1)*q+1:(ii+n-1)*q);
end
% for ii = 1:p
%     for jj = 1:q
%         H(ii:p:end,jj:q:end) = hankel(x(ii,jj:q:jj+q*(m-1)),x(ii,jj+q*(m-1):q:end));
%     end
% end


