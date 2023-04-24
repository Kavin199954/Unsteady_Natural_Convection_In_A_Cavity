%Name - Kavindu Jayasankha
%Roll No - 20BME116
%Topic - unsteady natural convection in cavity

clear all
clc

% Define problem parameters
n = 51;% % Number of grid points

x = linspace(0, 1, n);% X-coordinates of the grid
y = linspace(0, 1, n);% Y-coordinates of the grid 
dt = 0.001;           % time step for the simulation
endTime = 1;          % end time for the simulation

WriteInterval = 10; % Timesteps
h = x(2) - x(1); %grid spacing 
nu = 0.1; % Kinematic viscosity
alpha = 0.1;% thermal diffusivity
g = 9.81;  % Gravitational acceleration
T0 = 100; % temperature at the right 
T1 = 500;% temperature at the left
beta = 1 / T0;% inverse of the reference temperature

[xx, yy] = meshgrid(x, y);
vort = zeros(n, n);% matrix containing vorticity values at each grid point
psi = zeros(n, n);% streamfunction values at each grid point
temp = T0 * ones(n, n);% % Temperature field
temp(:, 1) = T1;

t = 0;
iter = 0;
run = 1;


% Add these lines before the main loop

u_history = [];
temp_history = [];

% Create a new VideoWriter object
writerObj = VideoWriter('output.mp4', 'MPEG-4');
writerObj.FrameRate = 30; % set the frame rate of the video

open(writerObj); % open the video writer object

while t < endTime
    vortold = vort;
    psiold = psi;
    tempold = temp;
    iter = iter + 1;
    t = t + dt;
    % Update vorticity
    for i = 2:n-1
        for j = 2:n-1
            u = (psi(i, j + 1) - psi(i, j - 1)) / (2 * h);
            v = (psi(i + 1, j) - psi(i - 1, j)) / (2 * h);
            vort(i, j) = vort(i, j) + dt * (-(psi(i, j + 1) - psi(i, j - 1)) * (vort(i + 1, j) - vort(i - 1, j)) / (4 * h ^ 2) ...
                + (psi(i + 1, j) - psi(i - 1, j)) * (vort(i, j + 1) - vort(i, j - 1)) / (4 * h ^ 2) ...
                + nu * (vort(i + 1, j) + vort(i, j + 1) + vort(i - 1, j) + vort(i, j - 1) - 4 * vort(i, j)) / h ^ 2 ...
                - beta * g * (temp(i, j+1) - temp(i, j-1)) / (2 * h));
        end
    end
    % Solve for streamfunction
    for i = 2:n-1
        for j = 2:n-1
            psi(i, j) = 0.25 * (h ^ 2 * vort(i, j) + psi(i + 1, j) + psi(i - 1, j) + psi(i, j + 1) + psi(i, j - 1));
        end
    end
    % Update temperature
    for i = 2:n-1
        for j = 2:n-1
            u = (psi(i, j + 1) - psi(i, j - 1)) / (2 * h);
            v = (psi(i + 1, j) - psi(i - 1, j)) / (2 * h);
            temp(i, j) = temp(i, j) + dt * (alpha * (temp(i + 1, j) + temp(i - 1, j) + temp(i, j + 1) + temp(i, j - 1) - 4 * temp(i, j)) / h ^ 2 ...
                - u * (temp(i, j + 1) - temp(i, j - 1)) / (2 * h) ...
                - v * (temp(i + 1, j) - temp(i - 1, j)) / (2 * h));
        end
    end
    % Impose temperature boundary conditions
    temp(:, 1) = T1; % Left
    temp(:, end) = T0; % Right
    temp(1, :) = temp(2, :); % Bottom
    temp(end, :) = temp(end - 1, :); % Top
    
    % Update velocities
    for i = 1:n
        for j = 1:n
            if j ~= 1 && j ~= n && i ~= 1 && i ~= n
                u(i, j) = (psi(i, j + 1) - psi(i, j - 1)) / (2 * h);
                v(i, j) = (psi(i + 1, j) - psi(i - 1, j)) / (2 * h);
            else
                u(i, j) = 0;
                v(i, j) = 0;
            end
        end
    end
    u_history = [u_history, u(:)];
    temp_history = [temp_history, temp(:)];
    if mod(iter, WriteInterval) == 0
        fprintf('Calculating ... time = %f \n', t)
    end
   % Temperature contour plot
   subplot(1,2,1)
   contourf(y, x, temp', 40, 'LineColor', 'none')
   axis equal
   axis([0 1 0 1])
   xlabel('x')
   ylabel('y')
   title('Temperature')
   colorbar  
    
   % Streamline plot for velocity
   subplot(1,2,2)
   quiver(x,y,u,v)
   axis equal
   axis([0 1 0 1])
   xlabel('x')
   ylabel('y')
   title('velocity')
    
    % Write the current frame to the video
    frame = getframe(gcf); % capture the current figure as a frame
    writeVideo(writerObj, frame); % write the frame to the video
end
close(writerObj); % close the video writer object

figure('name', 'Results')

% Streamline plot for velocity
subplot(2, 2, 1)
streamslice(y, x, u', v')
axis equal
axis([0 1 0 1])
xlabel('x')
ylabel('y')
title('Streamlines')

% Temperature contour plot
subplot(2, 2, 2)
contourf(y, x, temp', 40, 'LineColor', 'none')
axis equal
axis([0 1 0 1])
xlabel('x')
ylabel('y')
title('Temperature')
colorbar

% XY plot for velocity
subplot(2, 2, 3)
plot(linspace(0, t, size(u_history, 2)), max(u_history, [], 1))
xlabel('Time')
ylabel('Max velocity')
title('Velocity vs Time')

%XY plot for temperature
subplot(2, 2, 4)
plot(linspace(0, t, size(temp_history, 2)), max(temp_history, [], 1))
xlabel('Time')
ylabel('Max temperature')
title('Temperature vs Time')
