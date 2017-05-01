function [LUTpara] = create_LUTparam(turbi,turbines,wakes,LUTpara)
% This function calculates the parameters for the LUT.

    locTurbi      = turbines(turbi).LocWF;             % Determine location of current turbine
    largestImpact = turbines(turbi).turbLargestImpact; % Determine index of the turbine with largest impact on current turbine. 
    
    if isempty(largestImpact) == 0
        index = wakes(largestImpact).centerLine(1,:)==locTurbi(1);
        LUTpara.C2C = abs(locTurbi(2)-wakes(largestImpact).centerLine(2,index)); % Center to center distance
        LUTpara.Dw  = wakes(largestImpact).diameters(index,3);                   % Diameter of mixing zone
    else
        LUTpara.C2C = 230; % No overlap results in maximum C2C
        LUTpara.Dw  = 179; % No overlap results in minimum Dw
    end
    
    LUTpara.yaw  = turbines(turbi).YawIF;     % Current yaw angles
    LUTpara.Ueff = turbines(turbi).windSpeed; % Effective inflow speed
end