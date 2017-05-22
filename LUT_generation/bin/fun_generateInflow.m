function [] = fun_generateInflow(inputData,filename)

% Import settings
T                   = inputData.T;        
dt                  = inputData.dt;   
plotProfile         = inputData.plotProfile;
HH                  = inputData.HH;
yWake               = inputData.yWake;
zWake               = inputData.zWake;
yaw                 = inputData.yaw;
U_fs                = inputData.U_fs;     
Gaussian_A          = inputData.Gaussian_A;
Gaussian_omegay     = inputData.Gaussian_omegay;
Gaussian_omegaz     = inputData.Gaussian_omegaz;
Drotor              = inputData.Drotor; 
Dwake               = inputData.Dwake;
shear_const         = inputData.shear_const;
doWindShear         = inputData.doWindShear;


%% Parameters Formules Gebraad

ai   = 1/3;                 % Axial Induction Factor        [-]
k_e  = 0.065;               % Gebraad                       [-]
m_e  = [-0.5 0.22 1];       % Gebraad                       [-]
a_U  = 5;                   % Gebraad                       [-]
b_U  = 1.66;                % Gebraad                       [-]
M_U  = [0.5 1 5.5];         % Gebraad                       [-]
X    = 0;                   % Turbine location which creates the downstream flow [m]
q1   = 1;                   % Q value innermost wakezone
q3   = 3;                   % Q value outermost wakezone
yawt = 0;                   % Yaw upstream turbine          [deg]            

%% Static settings: vertical grid
y     = inputData.y;   % lateral dimension (NOTE: MUST BE POSITIVE TO NEGATIVE).
z     = inputData.z;   % vertical dimension
Ny    = length(y);     % Number of grid points y-
Nz    = length(z);     % Number of grid points z-
[Y,Z] = ndgrid(y,z);   % 2D grid points
time      = [dt:dt:T];     % Time vector [s]
x         = U_fs*time;     % longitudinal dimension [m]
Nx        = length(x);     % Number of grid points x-

%% Wake Zone Formula's Gebraad 
    
xwake = X + (Dwake - Drotor)/(2*k_e*m_e(q3));

%Get the value of Ueff coupled to the length in x-direction
[m_u,c,Ueff] = deal(zeros(1,3));    % initialize vectors
for q = 1:3
    m_u(q)  = M_U(q)/((cosd(a_U+b_U*yawt)));
    c(q)    = ((Drotor/(Drotor+2*k_e*m_u(q)*(xwake-X)))^2);
    Ueff(q) = U_fs*(1-2*ai*c(q));
 end
    
%%  Gaussian shape
Gaussian_omegay = Dwake/4;
Gaussian_omegaz = Dwake/4;

wakeGrid = zeros(Ny,Nz); % Calculate wake deficit
for dyi = 1:Ny
    dy = y(dyi)-yWake;
    for dzi = 1:Nz
        dz = z(dzi)-zWake;
        wakeGrid(dyi,dzi) =  ((U_fs-Ueff(q1))* exp(-(  ((dy.^2)/(2*Gaussian_omegay^2) + (dz.^2)/(2*Gaussian_omegaz^2))  )));
        % https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
    end
end

% Calculate one slice of the windfield
u_waked = U_fs*ones(Ny,Nz)-wakeGrid; % maybe add round(..,N)?
%% calculate windfield with shear effects

if(doWindShear)
    % Calculate shear effects
    v = zeros(1,length(z)-1);   % initialize vector
    for i = 2:length(z)
        v(i) = U_fs/((HH/z(i))^shear_const); % Velocity distribution in z-direction due to shear effects
    end

    % apply wind shear to windfield
    for i = 1:length(z)
        u_waked(:,i) = u_waked(:,i)*(v(i)/U_fs);    % multiply u_waked with normalized velocity distribution
    end
end

% generate u_out, v_out, w_out
u_fs = U_fs;
% Copy and add turbulence to the slices
TI = 0.0; % Currently 0. We have to think about time sampling and TI...
[u_out,v_out,w_out] = deal(zeros(Nx,Ny,Nz));
for i = 1:Nx
    u_out(i,:,:) = u_waked+u_waked*(TI*randn);
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
end

% Save to external files for FAST usage (.wnd)
% --- filename needs to be extended according to added dimensions to LUT ---
%filename = ['inflowProfiles/' destinationFolder '/' inflowFilename(inputData)];
writebladed(filename,(u_out-u_fs)/u_fs,v_out,w_out,x,y,z,U_fs);
fid = fopen([filename, '.sum'], 'wt'); % Write .sum file
fprintf(fid, 'T\tCLOCKWISE\n');
fprintf(fid, '%0.0f\tHUB HEIGHT\n\n', HH);
fprintf(fid, '%0.3f\tUBAR\n', U_fs);
fprintf(fid, '%0.3f\tTI(u)\n', 100);
fprintf(fid, '%0.3f\tTI(v)\n', 100);
fprintf(fid, '%0.3f\tTI(w)\n\n', 100);
fprintf(fid, '0\tHEIGHT OFFSET');
fclose(fid);
end