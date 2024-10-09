% Enter input parameters
L = 1.5;   % Length of interval
c = 1;   % Wave speed
tf = 20;  % Final time 
nf = 20; % Number of Fourier modes to keep
nx = 200; % Number of meshpoints to keep
nt = 1000;  % Number of time points to print

mesh = linspace(0,L,nx);
time = linspace(0,tf,nt);

% Calculate Fourier coefficients of initial conditions
% and calculate omega_n
for n=1:nf
    func = f0(mesh,L).*sin(n*pi*mesh/L);
    a(n) = (2/L)*trapz(mesh,func);
    omega(n) = n*pi*c/L;
end

plot(mesh,f0(mesh,L));

for i=1:nt
    for j=1:nx
        sum = 0;
        for n=1:nf
            sum = sum + a(n)*sin(n*pi*mesh(j)/L)*cos(omega(n)*time(i));
        end
        sol(j,i) = sum;
    end
    plot(mesh,sol(:,i));
    xlim([0 1]);
    ylim([-2 2]);
    M(i) = getframe;
end

movie(M);

function y=f0(x,L)
%y = sin(3*pi*x/L);
%y = x.*(L-x);
y = sin(3*pi*x/L).^2;
end
