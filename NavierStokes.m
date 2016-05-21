function [F] = NavierStokes()
%Navier Stokes Solver for the channel

%parameters we can change
L  = .1;
p  = 2000;
nu = .01;
dt = .0001;
T  = 10;
N  = 143;  %N must be odd , number of gridpoints in y direction
V  = 10;

%derived values from parameters
h = 1/(N-1); %size interval
K      = (N+1)/2;
S      = 2*ceil(L/(2*h)) + 1;
topS   = K + (S-1)/2;
botS   = K - (S-1)/2;
rightS = K/2 + S - 1;
x_size = rightS + 6*S;
K      = K/2;

oldVort = zeros(N, x_size);
newVort = zeros(N, x_size);
stream = zeros(N, x_size);

F(p/20) = struct('cdata',[],'colormap',[]);

figure
ax = gca;
ax.NextPlot = 'replaceChildren';

%initialize stream function at boundaries
for i = 1:N
    stream(i,1) = 4*V*((1/2)*((i-1)*h)^2 - (1/3)*((i-1)*h)^3);
end

stream(N, :) = 2*V/3;
stream(botS:topS, K:rightS) = V/3;
stream(botS+1:topS-1, K+1:rightS-1) = 0;

%build sparse matrix for inversion 
A = sparse(N*x_size, N*x_size);
A(1:N, 1:N) = eye(N);
A((x_size-1)*N+1 : x_size*N, (x_size-1)*N+1:x_size*N) = eye(N);
%newcode
A((x_size-1)*N+1 : x_size*N, (x_size-1)*N+1-N:x_size*N-N) = -eye(N);

for i= N+1 : (x_size-1)*N

    if mod(i,N) == 0 || mod(i-1,N) == 0
        A(i,i) = 1;
        
    elseif mod(i,N) >= botS && mod(i,N) <= topS && i-mod(i,N) >= N*(K-1) && i-mod(i,N) <= N*(rightS-1)
        A(i,i) = 1;
        
    else
    A(i, i-N)     = -dt*nu/(2*h*h);
    A(i, i-1:i+1) = [-dt*nu/(2*h*h), 1+2*dt*nu/(h*h), -dt*nu/(2*h*h)];
    A(i, i+N)     = -dt*nu/(2*h*h);
    end
end

b = zeros(N*x_size, 1);

for m = 1:p
    m
%SOR
w = 1.9;
totDiff = 1;
maxIter = 0;
while totDiff > .0001 && maxIter < 10000;
    totDiff = 0;
    for j = 2 : x_size-1
        for i = 2 : N-1
            if j < K || j > rightS || i < botS || i > topS
                temp = stream(i,j);
                stream(i,j) = w*(stream(i,j-1) + stream(i, j+1) + stream(i-1, j) + stream(i+1, j) + h^2*oldVort(i,j))/4 + (1-w)*stream(i,j); 
                totDiff = totDiff + (temp - stream(i,j))^2;
            end
        end
    end
    for i = 2 : N-1
        stream(i, x_size) = stream(i, x_size-1);
    end
    maxIter = maxIter + 1;
end

%set Vort at the outter boundaries
newVort(1,:) = 2*(stream(1, :) - stream(2, :))/(h*h);
newVort(N,:) = 2*(stream(N, :) - stream(N-1, :))/(h*h);
newVort(2:N-1, 1) = (stream(1:N-2, 1) + stream(3:N, 1) - 2*stream(2:N-1, 1))/(h*h) + 2*(stream(2:N-1, 1) - stream(2:N-1, 2))/(h*h);
%part of old code
%newVort(2:N-1, x_size) = (stream(1:N-2, x_size) + stream(3:N, x_size) - 2*stream(2:N-1, x_size))/(h*h) + 2*(stream(2:N-1, x_size) - stream(2:N-1, x_size-1))/(h*h);
%new code
newVort(1:N, x_size) = 0;

%set Vort at the square boundaries
newVort(botS+1:topS-1, K) = 2*(stream(botS+1:topS-1, K) - stream(botS+1:topS-1, K-1))/(h*h);
newVort(botS+1:topS-1, rightS) = 2*(stream(botS+1:topS-1, rightS) - stream(botS+1:topS-1, rightS+1))/(h*h);
newVort(botS, K+1:rightS-1) = 2*(stream(botS, K+1:rightS-1) - stream(botS-1, K+1:rightS-1))/(h*h);
newVort(topS, K+1:rightS-1) = 2*(stream(topS, K+1:rightS-1) - stream(topS+1, K+1:rightS-1))/(h*h);

newVort(topS,K) = (newVort(topS,K+1) + newVort(topS-1,K))/2;
newVort(botS,K) = (newVort(botS,K+1) + newVort(botS+1, K))/2;
newVort(topS,rightS) = (newVort(topS,rightS-1) + newVort(topS-1,rightS))/2;
newVort(botS,rightS) = (newVort(botS,rightS-1) + newVort(botS+1, rightS))/2;

for j = 2:x_size-1
    for i = 2 : N-1
        if j < K || j > rightS || i < botS || i > topS
            newVort(i,j) = oldVort(i,j) - dt*((stream(i+1,j) - stream(i-1,j))*(oldVort(i,j+1) - oldVort(i,j-1) + (dt/h)*(oldVort(i,j+1) -2*oldVort(i,j) + oldVort(i,j-1))) - (stream(i,j+1) - stream(i,j-1))*(oldVort(i+1,j) - oldVort(i-1,j)))/(4*h*h) + dt*nu*(oldVort(i,j+1) + oldVort(i,j-1) + oldVort(i+1,j) + oldVort(i-1,j) - 4*oldVort(i,j))/(2*h*h);
        end
    end
end
oldVort = reshape(A\(newVort(:)), N, x_size);

if mod(m,20) == 0
contourf(oldVort, linspace(-1000, 1000, 40));
axis equal
colormap('jet')
colorbar;
drawnow
F(m/20) = getframe(gcf);
end

end

fig = figure;
movie(fig, F, 1)
end














