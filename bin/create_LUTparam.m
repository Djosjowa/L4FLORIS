function [LUTpara] = create_LUTparam(turbi,turbines,wakes,LUTpara)
% This function calculates the parameters for the LUT.
% TODO: min and max parameter values are now static, make them variable

    locTurbi      = turbines(turbi).LocWF;             % Determine location of current turbine
    largestImpact = turbines(turbi).turbLargestImpact; % Determine index of the turbine with largest impact on current turbine.
    min_Dwake = 180;
    max_yWake = 250;
    min_U_fs = 6;
    max_U_fs = 8;
    min_yaw = -30;
    max_yaw = 30;
    
    % Get yWake and Dwake information
    if isempty(largestImpact) == 0
        index = wakes(largestImpact).centerLine(1,:)==locTurbi(1);
        LUTpara.yWake = abs(locTurbi(2)-wakes(largestImpact).centerLine(2,index)); % Center to center distance
        LUTpara.Dwake  = wakes(largestImpact).diameters(index,3);                   % Diameter of mixing zone
    else
        LUTpara.Dwake = min_Dwake; % No overlap results in minimum Dw
        LUTpara.yWake = max_yWake; % No overlap results in maximum C2C
    end
    
    % if Dwake or yWake falls outside the LUT boundaries, set it to boundary values (almost freestream conditions)
    if LUTpara.Dwake < min_Dwake || abs(LUTpara.yWake) > max_yWake
        LUTpara.Dwake = min_Dwake; % No overlap results in minimum Dw
        LUTpara.yWake = max_yWake; % No overlap results in maximum C2C
    end
    
    % Get yaw and freestream velocity information
    LUTpara.yaw  = turbines(turbi).YawIF;     % Current yaw angles
    LUTpara.U_fs = turbines(turbi).windSpeed; % Effective inflow speed
    
    % if U_fs and yaw are not in range, set them to min or max value
    if LUTpara.yaw > max_yaw
        LUTpara.yaw = max_yaw;
    end
    if LUTpara.yaw < min_yaw
        LUTpara.yaw = min_yaw;
    end
    if LUTpara.U_fs > max_U_fs
        LUTpara.U_fs = max_U_fs;
    end
    if LUTpara.U_fs < min_U_fs
        LUTpara.U_fs = min_U_fs;
    end
    
end