% Enter input parameters
L = 1;   % Length of interval
H = 2; 
nf = 20; % Number of Fourier modes to keep
nx = 50; % Number of meshpoints to keep in x 
ny = 100; % Number of meshpoints to keep in y 

meshx = linspace(0,L,nx);
meshy = linspace(0,H,ny);

% Calculate Fourier coefficients of initial conditions
% and calculate tau_n
for n=1:nf
    sqlambda(n) = n*pi/H;
    func = f0(meshy,H).*sin(n*pi*meshy/H);
    b(n) = (2/H)*trapz(meshy,func)/sinh(sqlambda(n)*L);
end

 for j=1:ny
    for i=1:nx
        sum = 0;
        for n=1:nf
            sum = sum + b(n)*sinh(sqlambda(n)*meshx(i))*sin(sqlambda(n)*meshy(j));
        end
        sol(j,i) = sum;
    end
end

surf(meshx,meshy,sol)
xlabel('x')
ylabel('y')

function res=f0(y,L)
%res = sin(3*pi*y/L);
%res = y.*(L-y);
res = sin(3*pi*y/L).^2;
end
