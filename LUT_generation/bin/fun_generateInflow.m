function [] = fun_generateInflow(inputData,filename)

% Import settings
T                   = inputData.T;        
dt                  = inputData.dt;   
plotProfile         = inputData.plotProfile;
HH                  = inputData.HH;   
yWake               = inputData.yWake;
zWake               = inputData.zWake;
yaw                 = inputData.yaw;
u_fs                = inputData.u_fs;     
Gaussian_A          = inputData.Gaussian_A;
Gaussian_omegay     = inputData.Gaussian_omegay;
Gaussian_omegaz     = inputData.Gaussian_omegaz;

% Static settings: vertical grid
y     = inputData.y;   % lateral dimension (NOTE: MUST BE POSITIVE TO NEGATIVE).
y_rot = y.*cosd(-yaw); % used instead of y, when there is a yaw angle. Compensates for resolution loss when rotating windfield
z     = inputData.z;   % vertical dimension
Ny    = length(y);     % Number of grid points y-
Nz    = length(z);     % Number of grid points z-
[Y,Z] = ndgrid(y,z);   % 2D grid points

time      = [dt:dt:T];     % Time vector [s]
x         = u_fs*time;     % longitudinal dimension [m]
Nx        = length(x);     % Number of grid points x-

wakeGrid = zeros(Ny,Nz); % Calculate wake deficit
for dyi = 1:Ny
    dy = y_rot(dyi)-yWake;
    for dzi = 1:Nz
        dz = z(dzi)-zWake;
        wakeGrid(dyi,dzi) = Gaussian_A * exp(-(  ((dy.^2)/(2*Gaussian_omegay^2) + (dz.^2)/(2*Gaussian_omegaz^2))  ));
        % https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
    end;
end;

% Calculate one slice of the windfield
u_waked = u_fs*ones(Ny,Nz)-wakeGrid; % maybe add round(..,N)?
if yaw == 0
    u_mean = u_fs; v_mean = 0;  % u_mean is the same as free stream velocity
    % Copy and add turbulence to the slices
    TI = 0.0; % Currently 0. We have to think about time sampling and TI...
    [u_out,v_out,w_out] = deal(zeros(Nx,Ny,Nz));
    for i = 1:Nx
        u_out(i,:,:) = u_waked+u_waked*(TI*randn);
    end;
    
else % if yaw angle isn't zero, rotate the windfield
    % calculate u_mean and v_mean
    u_mean = u_fs*cosd(-yaw);
    v_mean = u_fs*sind(-yaw);
    
    % Rotate u_waked
    rotu_waked = u_waked.*cosd(-yaw);
    rotv_waked = u_waked.*sind(-yaw);
    
    % Copy  the slices
    % Turbulence is ignored here for now, can be added later
    [u_out,v_out,w_out] = deal(zeros(Nx,Ny,Nz));
    for i = 1:Nx
        u_out(i,:,:) = rotu_waked;
        v_out(i,:,:) = rotv_waked;
    end;
end

% Plotting wake profile
if plotProfile
    
    % Plot front profile
    figure(1); clf; 
%     sp1 = subplot(2,1,1);
    contourf(Y,Z,reshape(u_out(1,:,:),[size(u_waked)]));
    axis equal; xlabel('y (m)'); ylabel('z (m)'); title('Inflow profile (m/s)');
    colorbar; zlabel('Velocity in x direction (m/s)'); hold on;
    plot(0,zWake,'r+');
    drawnow;
    
%     sp2 = subplot(2,1,2); contourf(Y,Z,u_waked);
%     axis equal; xlabel('y (m)'); ylabel('z (m)'); title('Inflow profile original (m/s)');
%     colorbar; zlabel('Velocity in x direction (m/s)'); hold on; linkaxes
%     plot(0,zWake,'r+'); linkaxes([sp1 sp2],'x');
%     drawnow;
%     set([sp1 sp2],'clim',[3 8]);
end;

% Save to external files for FAST usage (.wnd)
% --- filename needs to be extended according to added dimensions to LUT ---
%filename = ['inflowProfiles/' destinationFolder '/' inflowFilename(inputData)];
if(yaw==0)
    writebladed(filename,(u_out-u_mean)/u_mean,v_out,w_out,x,y,z,u_fs);
else
    writebladed(filename,(u_out-u_mean)/u_mean,(v_out-v_mean)/v_mean,w_out,x,y,z,u_fs);
end
fid = fopen([filename, '.sum'], 'wt'); % Write .sum file
fprintf(fid, 'T\tCLOCKWISE\n');
fprintf(fid, '%0.0f\tHUB HEIGHT\n\n', HH);
fprintf(fid, '%0.3f\tUBAR\n', u_fs);
fprintf(fid, '%0.3f\tTI(u)\n', 100);
fprintf(fid, '%0.3f\tTI(v)\n', 100);
fprintf(fid, '%0.3f\tTI(w)\n\n', 100);
fprintf(fid, '0\tHEIGHT OFFSET');
fclose(fid);
end