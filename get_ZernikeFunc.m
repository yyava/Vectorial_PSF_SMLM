function W = get_ZernikeFunc(N,aberrations)

n = aberrations(:,1);
m = aberrations(:,2);
c = aberrations(:,3)*2*pi;

d = 2/N;
xl = -(1-d/2):d:(1-d/2);

[X,Y] = meshgrid(xl,xl);

[theta,r] = cart2pol(Y,X);

idx = (r<=1);

z = zeros(N);
W = zeros(N);
y = zernfun(n,m,r(idx),theta(idx),'norm');

for i = 1:length(c)
    z(idx) = y(:,i)*c(i);
    W = W+z;
end
