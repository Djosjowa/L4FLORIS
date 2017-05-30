function[a_opt,yaw_opt,J_Pws_opt,J_DEL_opt,J_sum_opt] = optimizeL4FLORIS(modelStruct,turbType,siteStruct,optimStruct,LUT,Pref,Pbandwidth,plotResults)
% Optimization parameters
optConst    = optimStruct.optConst;
iterations  = optimStruct.iterations;
yawmin      = optimStruct.minYaw;
yawmax      = optimStruct.maxYaw;
yawinit     = optimStruct.initYaw;
amin        = optimStruct.minA;
amax        = optimStruct.maxA;
ainit       = optimStruct.initA;
N           = size(siteStruct.LocIF,1); % Number of turbines
DELbaseline = mean(mean(mean(mean(LUT.table)))); % DEL values are scaled with this value in the cost function
Pref_plot   = zeros(iterations,1)';

% Calculate windspeed distribution in wind-aligned frame
windSpeed                = hypot(siteStruct.uInfIf,siteStruct.vInfIf); % Static Wind Speed [m/s]
windDirection            = atand(siteStruct.vInfIf/siteStruct.uInfIf); % Nominal wind direction
windInflowDistribution   = windDirection+optimStruct.windUncertainty;  % Uncertain wind directions
weightsInflowUncertainty = gaussianWindDistribution(windInflowDistribution,plotResults); % Weights for inflow

% Initialize empty GT-theory matrices
[J_Pws_opt,J_sum_opt] = deal(-1e10);
J_DEL_opt             = 1e10;
yaw                   = yawinit*ones(N,1);
yaw_opt               = yawinit*ones(iterations,N);
a                     = ainit*ones(N,1);
a_opt                 = ainit*ones(iterations,N); 

% Perform game-theoretic optimization
disp([datestr(rem(now,1)) ': Starting GT optimization using FLORIS. [Iterations: ' num2str(iterations) '. Calls to FLORIS: ' num2str(iterations*length(windInflowDistribution)) ']']); tic;
for k = 1:iterations  % k is the number of iterations
    if(~rem(k*100/iterations,10)); disp([datestr(rem(now,1)) ':  ' num2str(k*100/iterations) '% completed.']); end;
    
    % For k == 1 do a baseline run, otherwise randomize yaw angles
    if k > 1
        for i = 1:N                 % For each WT
            R1 = rand();            % Random value between [0 1]
            R2 = rand();
            E = 1-k/iterations;     % Sensitivity linearly related to iteration
            if R1 < E
                R3 = normrnd(0,0.2); % Perturb with random value
                a(i) = max(min(a_opt(i)+R3,amax),amin);
            else
                a(i) = a_opt(k-1,i);  
            end
            if R2 < E
                R4 = normrnd(0,35);
                yaw(i) = max(min(yaw_opt(i)+R4,yawmax),yawmin);
            else
                yaw(i) = yaw_opt(k-1,i);
            end
        end
    end
    
    for jj = 1:length(windInflowDistribution) % Calculate power and DEL for each wind direction
        windDir = windInflowDistribution(jj);
        siteStruct.uInfIf = windSpeed*cosd(windDir);
        siteStruct.vInfIf = windSpeed*sind(windDir);
        
        input.a   = a;
        input.yaw = yaw;
        
        [turbines,wakes] = run_floris(input,modelStruct,turbType,siteStruct);
        
        [P,DEL]  = deal(zeros(1,N));
        LUTparam = struct('Dwake',[],'U_fs',[],'yaw',num2cell(zeros(1,N)),'yWake',[]);
        
        for turbi = 1:N
            LUTparam(turbi) = create_LUTparam(turbi,turbines,wakes,LUTparam(turbi));
            P(turbi) = turbines(turbi).power;
            % -- Look up DEL values for flow field with C2C, Dw, Ueff
            DEL(turbi)= interpn(LUT.Dwake,LUT.U_fs,LUT.yaw,LUT.yWake,...
                                LUT.table,LUTparam(turbi).Dwake,LUTparam(turbi).U_fs,LUTparam(turbi).yaw,LUTparam(turbi).yWake); 
        end
        
        Ptot   = sum(P);
        DELtot = sum(DEL);
        
        Ptot_inflows(1,jj)   = Ptot;   % Store results for each wind direction
        DELtot_inflows(1,jj) = DELtot; % Store results for each wind direction
    end
    
    % Calculate collective results over entire wind rose
    sum_Ptot    = Ptot_inflows   * weightsInflowUncertainty;  % Inflow uncertainty-weighed generated power
    sum_DELtot  = DELtot_inflows * weightsInflowUncertainty;  % Inflow uncertainty-weighed turbine DEL values
    sum_PDELtot = optConst*((Pref-sum_Ptot)/Pbandwidth)^2 + (1-optConst)*sum_DELtot/DELbaseline; % Generate combined power and loads cost function
    
    if (sum_PDELtot < J_sum_opt | k == 1)
        a_opt(k,:)       = a;
        yaw_opt(k,:)     = yaw;
        J_Pws_opt(k)     = sum_Ptot;
        J_DEL_opt(k)     = sum_DELtot;
        J_sum_opt(k)     = sum_PDELtot;
        J_sum_sub_DEL(k) = (1-optConst)*sum_DELtot/DELbaseline;
        J_sum_sub_P(k)   = optConst*((Pref-sum_Ptot)/Pbandwidth)^2;
    else % if no improvements: keep optimal solution
        a_opt(k,:)       = a_opt(k-1,:);
        yaw_opt(k,:)     = yaw_opt(k-1,:);
        J_Pws_opt(k)     = J_Pws_opt(k-1);
        J_DEL_opt(k)     = J_DEL_opt(k-1);
        J_sum_opt(k)     = J_sum_opt(k-1);
        J_sum_sub_DEL(k) = J_sum_sub_DEL(k-1);
        J_sum_sub_P(k)   = J_sum_sub_P(k-1);
    end
    
    Pref_plot(k) = Pref;
    
    if k > 0.05*iterations
        J_sum_opt_95(k)     = J_sum_opt(k);
        J_sum_sub_DEL_95(k) = J_sum_sub_DEL(k);
        J_sum_sub_P_95(k)   = J_sum_sub_P(k);
    else
        J_sum_opt_95(k)     = NaN;
        J_sum_sub_DEL_95(k) = NaN;
        J_sum_sub_P_95(k)   = NaN;
    end
    
end
disp([datestr(rem(now,1)) ': Elapsed time is ' num2str(toc) ' seconds.']);

%% Plotting results
if plotResults
    disp([datestr(rem(now,1)) ': Plotting results...'])
        disp(a_opt(k,:));
        disp(yaw_opt(k,:));
    
    figure 
    % Cost function
    fig1 = subplot(2,2,1);
    hold on;
    plot(J_sum_opt,'Linewidth',2); grid on;
    plot(J_sum_sub_DEL);
    plot(J_sum_sub_P);
    title('Mixed optimization: Power & Load');
    ylabel('PL-score [-]'); xlabel('Iterations [-]');
    legend('PL-score','DEL part','Power part');
    
    fig2 = subplot(2,2,3);
    hold on
    plot(J_sum_opt_95,'Linewidth',2); grid on;
    plot(J_sum_sub_DEL_95);
    plot(J_sum_sub_P_95);
    title('Mixed optimization: Power & Load, last 95% of iterations');
    ylabel('PL-score [-]'); xlabel('Iterations [-]');
    legend('PL-score','DEL part','Power part');
    
    % Summed generated power
    fig3 = subplot(2,2,2);
    hold on
    plot(J_Pws_opt/1E6,'Linewidth',2);
    plot(Pref_plot/1E6,'--','Linewidth',1);
    grid on; title('Summed power');
    ylabel('Power [MW]'); xlabel('Iterations [-]');
    
    % Summed DEL values
    fig4 = subplot(2,2,4);
    plot(J_DEL_opt/1E6,'Linewidth',2)
    grid on; title('Summed DEL values')
    ylabel('DEL (10^6)'); xlabel('Iterations [-]')
    
    linkaxes([fig1 fig2 fig3 fig4],'x')
    
    % evaluate nominal wind direction, optimal settings, and plot output
    siteStruct.uInfIf      = windSpeed*cosd(windDirection);
    siteStruct.vInfIf      = windSpeed*sind(windDirection);
    plots.plotLayout = false ; % plot farm layout w.r.t. inertial and wind frame
    plots.plot2DFlow = true  ; % 2DflowFieldvisualisation in wind-aligned frame
    plots.plot3DFlow = false ; % 3DflowFieldvisualisation in wind-aligned frametimer.script = tic;
    run_floris(input,modelStruct,turbType,siteStruct,plots);
end
end