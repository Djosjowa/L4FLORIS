function [sum_PDELtot] = optimfunc(put,optimStruct,siteStruct,modelStruct,turbType,DEL_table,plotResults)

% Optimization parameters
N           = size(siteStruct.LocIF,1);             % Number of turbines
optConst    = optimStruct.optConst;
Pref        = optimStruct.Pref;
Pbandwidth  = optimStruct.Pbandwidth;
DELbaseline = mean(mean(mean(DEL_table.table)));    % DEL values are scaled with this value in the cost function
a     = put(1:9);
yaw   = put(10:18);
input.a = a;
input.yaw = yaw;
% Calculate windspeed distribution in wind-aligned frame
windSpeed                = hypot(siteStruct.uInfIf,siteStruct.vInfIf); % Static Wind Speed [m/s]
windDirection            = atand(siteStruct.vInfIf/siteStruct.uInfIf); % Nominal wind direction
windInflowDistribution   = windDirection+optimStruct.windUncertainty;  % Uncertain wind directions
weightsInflowUncertainty = gaussianWindDistribution(windInflowDistribution,plotResults); % Weights for inflow

for j = 1:length(windInflowDistribution) % Calculate power and DEL for each wind direction
    windDir = windInflowDistribution(j);
    siteStruct.uInfIf = windSpeed*cosd(windDir);
    siteStruct.vInfIf = windSpeed*sind(windDir);

    [turbines, wakes, wtRows] = run_floris(input,modelStruct,turbType,siteStruct);

    [P,DEL]  = deal(zeros(1,N));
    LUTparam = struct('C2C',[],'Dw',[],'yaw',num2cell(zeros(1,N)),'Ueff',[]);

    for turbi = 1:N
        LUTparam(turbi) = create_LUTparam(turbi,turbines,wakes,LUTparam(turbi));
        P(turbi) = turbines(turbi).power;
        % -- Look up DEL values for flow field with C2C, Dw, Ueff
        DEL(turbi)= interpn(DEL_table.C2C,DEL_table.Dw,DEL_table.Ueff,...
                            DEL_table.table,LUTparam(turbi).C2C,LUTparam(turbi).Dw,LUTparam(turbi).Ueff); 
        if isnan(DEL(turbi)) == 1   % This is a quick fix because values sometimes fall outside the LUT range
            DEL(turbi) = 0; 
        end
        % DEL(turbi) = 1; % Placeholder
    end

    Ptot   = sum(P);
    DELtot = sum(DEL);

    Ptot_inflows(1,j)   = Ptot;   % Store results for each wind direction
    DELtot_inflows(1,j) = DELtot; % Store results for each wind direction
end

% Costfunction
sum_Ptot    = Ptot_inflows   * weightsInflowUncertainty;  % Inflow uncertainty-weighed generated power
sum_DELtot  = DELtot_inflows * weightsInflowUncertainty;  % Inflow uncertainty-weighed turbine DEL values
sum_PDELtot = optConst*((Pref-sum_Ptot)/Pbandwidth)^2 + (1-optConst)*sum_DELtot/DELbaseline; % Generate combined power and loads cost function
