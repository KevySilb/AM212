% Enter input parameters
L = 1;   % Length of interval
D = 1;   % Diffusion coefficient 
tf = 1;  % Final time 
nf = 20; % Number of Fourier modes to keep
nx = 200; % Number of meshpoints to keep
nt = 100;  % Number of time points to print

mesh = linspace(0,L,nx);
time = linspace(0,tf,nt);

% Calculate Fourier coefficients of initial conditions
% and calculate tau_n
a0 = (1/L) * trapz(mesh,f0(mesh,L));
for n=1:nf
    func = f0(mesh,L).*cos(n*pi*mesh/L);
    a(n) = (2/L)*trapz(mesh,func);
    lambda(n) = D*n^2*pi^2/L^2;
end

plot(mesh,f0(mesh,L));
hold on;

for i=1:nt
    for j=1:nx
        sum = a0;
        for n=1:nf
            sum = sum + a(n)*cos(n*pi*mesh(j)/L)*exp(-lambda(n)*time(i));
        end
        sol(j,i) = sum;
    end
plot(mesh,sol(:,i));
end

function y=f0(x,L)
y = sin(3*pi*x/L);
%y = x.*(L-x);
%y = sin(3*pi*x/L).^2;
end
