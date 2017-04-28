function[yaw_opt,J_Pws_opt,J_DEL_opt,J_sum_opt] = optimizeL4FLORIS(modelStruct,turbType,siteStruct,optimStruct,plotResults )
% Optimization parameters
optConst   = optimStruct.optConst;
iterations = optimStruct.iterations;
input.a    = optimStruct.axInd;
yawmin     = optimStruct.minYaw;
yawmax     = optimStruct.maxYaw;
N          = size(siteStruct.LocIF,1); % Number of turbines

% Calculate windspeed distribution in wind-aligned frame
windSpeed                = hypot(siteStruct.uInfIf,siteStruct.vInfIf); % Static Wind Speed [m/s]
windDirection            = atand(siteStruct.vInfIf/siteStruct.uInfIf); % Nominal wind direction
windInflowDistribution   = windDirection+optimStruct.windUncertainty;  % Uncertain wind directions
weightsInflowUncertainty = gaussianWindDistribution(windInflowDistribution, plotResults); % Weights for inflow

% Initialize empty GT-theory matrices
[J_Pws_opt,J_sum_opt] = deal(-1e10);
J_DEL_opt             = 1e10;
yaw                   = zeros(N,1);           
yaw_opt               = zeros(iterations,N); 

% Perform game-theoretic optimization
disp([datestr(rem(now,1)) ': Starting GT optimization using FLORIS. [Iterations: ' num2str(iterations) '. Calls to FLORIS: ' num2str(iterations*length(windInflowDistribution)) ']']); tic;
for k = 1:iterations  % k is the number of iterations
    if(~rem(k*100/iterations,10)); disp([datestr(rem(now,1)) ':  ' num2str(k*100/iterations) '% completed.']); end;
    
    % For k == 1 do a baseline run, otherwise randomize yaw angles
    if k > 1
        for i = 1:N                 % For each WT
            R1 = rand();            % Random value between [0 1]
            E = 1-k/iterations;     % Sensitivity linearly related to iteration
            if R1 < E
                R2 = normrnd(0,35); % Perturb with random value [0 10] degrees
                yaw(i) = max(min(yaw_opt(i)+R2,yawmax),yawmin);
            else
                yaw(i) = yaw_opt(k-1,i);   % Keep old yaw setting [degrees]
            end;
        end;
    end;
    
    for jj = 1:length(windInflowDistribution) % Calculate power and DEL for each wind direction
        windDir = windInflowDistribution(jj);
        siteStruct.uInfIf = windSpeed*cosd(windDir);
        siteStruct.vInfIf = windSpeed*sind(windDir);
        
        input.yaw  = yaw;
        [turbines, wakes, wtRows] = run_floris(input,modelStruct,turbType,siteStruct);
        
        
        [P,DEL] = deal(zeros(1,N));
        for turbi = 1:N
            % Calculate parameters for LUT
            locTurbi = turbines(turbi).LocWF;
            largestImpact = turbines(turbi).turbLargestImpact;
            if isempty(largestImpact) == 0
                index = wakes(largestImpact).centerLine(1,:)==locTurbi(1);
                LUTparam(turbi).C2C = abs(locTurbi(2)-wakes(largestImpact).centerLine(2,index));
                LUTparam(turbi).Dw = wakes(largestImpact).diameters(index,3);
            end
            LUTparam(turbi).yaw = turbines(turbi).YawIF;
            LUTparam(turbi).Ueff = turbines(turbi).windSpeed;
            P(turbi)  = turbines(turbi).power;
            % -- Look up DEL values for flow field with value1, value2, value3 --
            % DEL(turbi)= interpn(DEL_table.param1,DEL_table.param2,DEL_table.param3,...
            %                    DEL_table.table,value1,value2,value3); 
            DEL(turbi) = 1; % Placeholder
        end;
        Ptot   = sum(P);
        DELtot = sum(DEL);
        
        if k==1 % Set baseline values for k == 1
            Pbaseline   = Ptot;
            DELbaseline = DELtot;
        end;
        
        Ptot_inflows(1,jj)   = Ptot;   % Store results for each wind direction
        DELtot_inflows(1,jj) = DELtot; % Store results for each wind direction
    end;
    
    % Calculate collective results over entire wind rose
    sum_Ptot    = Ptot_inflows   * weightsInflowUncertainty;  % Inflow uncertainty-weighed generated power
    sum_DELtot  = DELtot_inflows * weightsInflowUncertainty;  % Inflow uncertainty-weighed turbine DEL values
    sum_PDELtot = optConst*sum_Ptot/Pbaseline - (1-optConst)*sum_DELtot/DELbaseline; % Generate combined power and loads cost function
    
    if (sum_PDELtot >= J_sum_opt | k == 1)
        yaw_opt(k,:) = yaw;
        J_Pws_opt(k) = sum_Ptot;
        J_DEL_opt(k) = sum_DELtot;
        J_sum_opt(k) = sum_PDELtot;
    else % if no improvements: keep optimal solution
        yaw_opt(k,:) = yaw_opt(k-1,:);
        J_Pws_opt(k) = J_Pws_opt(k-1);
        J_DEL_opt(k) = J_DEL_opt(k-1);
        J_sum_opt(k) = J_sum_opt(k-1);
    end;
end;
disp([datestr(rem(now,1)) ': Elapsed time is ' num2str(toc) ' seconds.']);
if plotResults % Plot results
    disp([datestr(rem(now,1)) ': Plotting results...'])
    figure % Cost function
    plot(J_sum_opt,'Linewidth',2); grid on;
    title('Mixed optimization: Power & Load');
    ylabel('PL-score [-]'); xlabel('Iterations [-]');
    
    figure % Summed generated power
    plot(J_Pws_opt/1E6,'Linewidth',2);
    grid on; title('Summed power');
    ylabel('Power [MW]'); xlabel('Iterations [-]');
    
    figure % Summed DEL values
    plot(J_DEL_opt/1E6,'Linewidth',2)
    grid on; title('Summed DEL values')
    ylabel('DEL (10^6)'); xlabel('Iterations [-]')
    
    % evaluate nominal wind direction, optimal settings, and plot output
    siteStruct.uInfIf      = windSpeed*cosd(windDirection);
    siteStruct.vInfIf      = windSpeed*sind(windDirection);
    plots.plotLayout = false ; % plot farm layout w.r.t. inertial and wind frame
    plots.plot2DFlow = true  ; % 2DflowFieldvisualisation in wind-aligned frame
    plots.plot3DFlow = false ; % 3DflowFieldvisualisation in wind-aligned frametimer.script = tic;
    run_floris(input,modelStruct,turbType,siteStruct,plots);
end;
end