function [] = fun_generateInflow(inputData,filename)

% Import settings
T                   = inputData.T;        
dt                  = inputData.dt;   
plotProfile         = inputData.plotProfile;    
HH                  = inputData.HH;   
yWake               = inputData.yWake;
zWake               = inputData.zWake;
yaw                 = inputData.yaw;
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

% Calculate one slice of the windfield
u_waked = u_mean*ones(Ny,Nz)-wakeGrid; % maybe add round(..,N)?

if yaw == 0
    % Copy and add turbulence to the slices
    TI = 0.0; % Currently 0. We have to think about time sampling and TI...
    [u_out,v_out,w_out] = deal(zeros(Nx,Ny,Nz));
    for i = 1:Nx
        u_out(i,:,:) = u_waked+u_waked*(TI*randn);
    end;
    
else % if yaw angle isn't zero, rotate the windfield
    % Rotate u_waked and the Y and Z matrix
    % initialization of vectors for the rotated windspeed vectors
    rotu_waked = zeros(size(u_waked));
    rotv_waked = zeros(size(u_waked));

    % for each point in the yz-grid the location and the windspeed vectors
    % are rotated
    for iy = 1:length(y)
        for iz = 1:length(z)
            rotXYZ = rotz(yaw)*[0 ; Y(iy,iz) ; Z(iy,iz)];
            rotuvw = rotz(yaw)*[u_waked(iy,iz) ; 0 ; 0];
            Y(iy,iz) = rotXYZ(2);
            Z(iy,iz) = rotXYZ(3);
            rotu_waked(iy,iz) = rotuvw(1);
            rotv_waked(iy,iz) = rotuvw(2);
        end
    end
    
    % Copy  the slices
    % Turbulence is ignored here for now, can be added later
    [u_out,v_out,w_out] = deal(zeros(Nx,Ny,Nz));
    for i = 1:Nx
        u_out(i,:,:) = rotu_waked;
        v_out(i,:,:) = rotv_waked;
    end;
end

% save('plot_data','u_mean','wakeGrid','u_waked','u_out','v_out','w_out','Y','Z','Yx','X','x','y','z');

% Plotting wake profiles
if plotProfile
    % Plot front profile
    figure(1); clf; 
    subplot(2,1,1); contourf(Y,Z,reshape(u_out(1,:,:),[size(u_waked)]));
    axis equal; xlabel('y (m)'); ylabel('z (m)'); title('Inflow profile front view (m/s)');
    colorbar; zlabel('Velocity in x direction (m/s)'); hold on;
    plot(0,zWake,'r+');
    drawnow;
    
    % plot top profile
    % TODO: Yx isn't updated yet when windfield is rotated, so the y-axis of this
    % plot is wrong then
    subplot(2,1,2); contourf(Yx,X,u_out(:,:,z==zWake)');
    xlabel('y (m)'); ylabel('x (m)'); title('Inflow profile topview (m/s)');
    colorbar; zlabel('Velocity in x direction (m/s)'); hold on;
    drawnow;
end;

% Save to external files for FAST usage (.wnd)
% --- filename needs to be extended according to added dimensions to LUT ---
%filename = ['inflowProfiles/' destinationFolder '/' inflowFilename(inputData)];
writebladed(filename,(u_out-u_mean)/u_mean,v_out,w_out,x,Y(:,1)',z,u_mean); % Y(:,1)' instead of y is given to still be correct when windfield is rotated
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