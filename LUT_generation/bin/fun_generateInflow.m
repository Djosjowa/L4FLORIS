function [] = fun_generateInflow(inputData,filename)

% Import settings
T                   = inputData.T;        
dt                  = inputData.dt;   
plotProfile         = inputData.plotProfile;    
HH                  = inputData.HH;   
yWake               = inputData.yWake;
zWake               = inputData.zWake;           
u_mean              = inputData.u_mean;     
Gaussian_A          = inputData.Gaussian_A;
Gaussian_omegay     = inputData.Gaussian_omegay;
Gaussian_omegaz     = inputData.Gaussian_omegaz;

% Static settings: vertical grid
y     = inputData.y;   % lateral dimension (NOTE: MUST BE POSITIVE TO NEGATIVE)
z     = inputData.z;   % vertical dimension
Ny    = length(y);     % Number of grid points y-
Nz    = length(z);     % Number of grid points z-
[Y,Z] = ndgrid(y,z);   % 2D grid points

time      = [dt:dt:T];     % Time vector [s]
x         = u_mean*time;   % longitudinal dimension [m]
Nx        = length(x);     % Number of grid points x-
[Yx,X]     = ndgrid(y,x);    % 2D grid points

wakeGrid = zeros(Ny,Nz); % Calculate wake deficit
for dyi = 1:Ny
    dy = y(dyi)-yWake;
    for dzi = 1:Nz
        dz = z(dzi)-zWake;
        wakeGrid(dyi,dzi) = Gaussian_A * exp(-(  ((dy.^2)/(2*Gaussian_omegay^2) + (dz.^2)/(2*Gaussian_omegaz^2))  ));
        % https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
    end;
end;
% Plot inflow profile (m/s)
u_waked = u_mean*ones(Ny,Nz)-wakeGrid; % maybe add round(..,N)?

if plotProfile
    figure(1); clf; contourf(Y,Z,u_waked);
    axis equal; xlabel('y (m)'); ylabel('z (m)'); title('Inflow profile (m/s)');
    colorbar; zlabel('Velocity (m/s)'); hold on;
    plot(0,zWake,'r+');
    drawnow;
end;

% Copy and add turbulence to the slices
TI = 0.01; % Currently 0. We have to think about time sampling and TI...
[u_out,v_out,w_out] = deal(zeros(Nx,Ny,Nz));
for i = 1:Nx
    u_out(i,:,:) = u_waked+u_waked*(TI*randn);
end;

% Plot topview
if plotProfile
    % generating topview matrix
    
    % plotting
    figure(2); clf; contourf(Yx,X,u_out(:,:,find(z==zWake))');
    axis equal; xlabel('y (m)'); ylabel('x (m)'); title('Inflow profile topview (m/s)');
    colorbar; zlabel('Velocity (m/s)'); hold on;
    drawnow;
end;

% Save to external files for FAST usage (.wnd)
% --- filename needs to be extended according to added dimensions to LUT ---
%filename = ['inflowProfiles/' destinationFolder '/' inflowFilename(inputData)];
writebladed(filename,(u_out-u_mean)/u_mean,v_out,w_out,x,y,z,u_mean);
fid = fopen([filename, '.sum'], 'wt'); % Write .sum file
fprintf(fid, 'T\tCLOCKWISE\n');
fprintf(fid, '%0.0f\tHUB HEIGHT\n\n', HH);
fprintf(fid, '%0.3f\tUBAR\n', u_mean);
fprintf(fid, '%0.3f\tTI(u)\n', 100);
fprintf(fid, '%0.3f\tTI(v)\n', 100);
fprintf(fid, '%0.3f\tTI(w)\n\n', 100);
fprintf(fid, '0\tHEIGHT OFFSET');
fclose(fid);
end