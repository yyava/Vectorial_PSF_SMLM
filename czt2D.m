function Out = czt2D(In,M,w,a)
% This function evaluates the 2D FT via the czt-algorithm to the first 2
% dimensions. The other dimension will remain the same as the Input.
% Input: 
%   In size:(N,N,...)
% Output:
%   Out size(M,M,...)

s = size(In);
N = s(1);

K = prod(s(3:end));

In = reshape(In,[N,N,K]);

InterTerm = permute(czt(In,M,w,a),[2,1,3]);

Out = permute(czt(InterTerm,M,w,a),[2,1,3]);

Out = reshape(Out,[M,M,s(3:end)]);

end