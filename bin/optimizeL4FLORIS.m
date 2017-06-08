function[output] = optimizeL4FLORIS(modelStruct,turbType,siteStruct,optS,LUT,plotResults)
%% Initialization
% Optimization parameters
N           = size(siteStruct.LocIF,1);  % Number of turbines
it          = optS.iterations;           % Number of iterations
DELbaseline = mean2(LUT.table);          % Baseline DEL value for cost function

% Calculate windspeed distribution in wind-aligned frame
windSpeed                = hypot(siteStruct.uInfIf,siteStruct.vInfIf);                      % Static Wind Speed [m/s]
windDirection            = atand(siteStruct.vInfIf/siteStruct.uInfIf);                      % Nominal wind direction
windInflowDistribution   = windDirection+optS.windUncertainty;                              % Uncertain wind directions
weightsInflowUncertainty = gaussianWindDistribution(windInflowDistribution,plotResults);    % Weights for inflow

% Initialize empty GT-theory matrices
[J_Pws_opt,J_DEL_opt,J_sum_sub_DEL,J_sum_sub_P,...
    J_sum_opt_95,J_sum_sub_DEL_95,J_sum_sub_P_95]   = deal(zeros(1,it));
[yaw_opt,yaw_tries,a_opt,a_tries]                   = deal(zeros(it,N));
[Ptot_inflows,DELtot_inflows]                       = deal(zeros(1,length(windInflowDistribution)));

J_sum_opt = -1e10;
yaw       = optS.initYaw*ones(N,1);
a         = optS.initA*ones(N,1);
Pref_plot = zeros(it,1)';

%% Game-theoretic optimization
disp([datestr(rem(now,1)) ': Starting GT optimization using FLORIS. [Iterations: ' num2str(it) '. Calls to FLORIS: ' num2str(it*length(windInflowDistribution)) ']']); tic;
for k = 1:it  % k is the iteration number of the GT optimization
    if(~rem(k*100/it,10)); disp([datestr(rem(now,1)) ':  ' num2str(k*100/it) '% completed.']); end
    
    % For k == 1 a baseline run is executed, for k > 1 randomize yaw angles
    if k > 1
        for i = 1:N                  % For each WT
            R1 = rand();             % Uniform distributed random value between [0 1]
            R2 = rand();
            E = 1-k/it;              % Sensitivity linearly related to current iteration
            if R1 < E                % Take a step if R1 < the sensitivity
                R3 = normrnd(0,0.1); % Perturb with normally distributed random value
                a(i) = max(min(a_opt(k-1,i)+R3,optS.maxA),optS.minA); % Perturb axial induction factor, with stepsize R3, from the latest baseline value
            else
                a(i) = a_opt(k-1,i); % Save latest baseline value
            end
            if R2 < E
                R4 = normrnd(0,15);
                yaw(i) = max(min(yaw_opt(k-1,i)+R4,optS.maxYaw),optS.minYaw); % Perturb yaw angle, with stepsize R4, from the latest baseline value
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
        
        [turbines,wakes] = run_floris(input,modelStruct,turbType,siteStruct);   % Execute FLORIS model
        
        [P,DEL]  = deal(zeros(1,N));
        LUTparam = struct('Dwake',[],'U_fs',[],'yaw',num2cell(zeros(1,N)),'yWake',[]);
        
        for turbi = 1:N                                                                 % For each turbine
            P(turbi)        = turbines(turbi).power;                                    % Power as calculated by the FLORIS model [N]
            
            LUTparam(turbi) = create_LUTparam(turbi,turbines,wakes,LUTparam(turbi));    % Extract LUT parameters from the FLORIS model
            DEL(turbi)      = interpn(LUT.Dwake,LUT.U_fs,LUT.yaw,LUT.yWake,...          % Linear interpolation in LookUp Table
                                LUT.table,LUTparam(turbi).Dwake,LUTparam(turbi).U_fs,LUTparam(turbi).yaw,LUTparam(turbi).yWake); 
        end
        
        Ptot   = sum(P);
        DELtot = sum(DEL);
        
        Ptot_inflows(jj)   = Ptot;   % Store results for each wind direction
        DELtot_inflows(jj) = DELtot; % Store results for each wind direction
    end
    
    % Calculate collective results over entire wind rose
    sum_Ptot    = Ptot_inflows   * weightsInflowUncertainty; % Inflow uncertainty-weighed generated power
    sum_DELtot  = DELtot_inflows * weightsInflowUncertainty; % Inflow uncertainty-weighed turbine DEL values
    sum_PDELtot = optS.optConst*((optS.Pref-sum_Ptot)/optS.Pbandwidth)^2 +...
        (1-optS.optConst)*sum_DELtot/DELbaseline;            % Generate combined power and loads cost function
    
    if (k == 1 | sum_PDELtot < J_sum_opt)
        a_opt(k,:)       = a;
        yaw_opt(k,:)     = yaw;
        J_Pws_opt(k)     = sum_Ptot;
        J_DEL_opt(k)     = sum_DELtot;
        J_sum_opt(k)     = sum_PDELtot;
        J_sum_sub_DEL(k) = (1-optS.optConst)*sum_DELtot/DELbaseline;
        J_sum_sub_P(k)   = optS.optConst*((optS.Pref-sum_Ptot)/optS.Pbandwidth)^2;
    else % if no improvements: keep optimal solution
        a_opt(k,:)       = a_opt(k-1,:);
        yaw_opt(k,:)     = yaw_opt(k-1,:);
        J_Pws_opt(k)     = J_Pws_opt(k-1);
        J_DEL_opt(k)     = J_DEL_opt(k-1);
        J_sum_opt(k)     = J_sum_opt(k-1);
        J_sum_sub_DEL(k) = J_sum_sub_DEL(k-1);
        J_sum_sub_P(k)   = J_sum_sub_P(k-1);
    end
    
    a_tries(k,:) = a;
    yaw_tries(k,:) = yaw;
    Pref_plot(k) = optS.Pref;
    
    if k > 0.05*it
        J_sum_opt_95(k)     = J_sum_opt(k);
        J_sum_sub_DEL_95(k) = J_sum_sub_DEL(k);
        J_sum_sub_P_95(k)   = J_sum_sub_P(k);
    else
        J_sum_opt_95(k)     = NaN;
        J_sum_sub_DEL_95(k) = NaN;
        J_sum_sub_P_95(k)   = NaN;
    end
    
end

%% Plotting results
if plotResults
    disp(' ')
    disp('Optimized axial induction factors:')
    disp(a_opt(k,:))
    disp('Optimized yaw angles:')
    disp(yaw_opt(k,:))
    disp('Plotting results...')
    
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
disp([datestr(rem(now,1)) ': Elapsed time is ' num2str(toc) ' seconds.']);

% Send desired output to workspace
output.a_opt     = a_opt;
output.a_tries   = a_tries;
output.yaw_opt   = yaw_opt;
output.yaw_tries = yaw_tries;
output.J_Pws_opt = J_Pws_opt;
output.J_DEL_opt = J_DEL_opt;
output.J_sum_opt = J_sum_opt;

end