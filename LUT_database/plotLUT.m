%% Settings
clear; close all;
LUTname = '4D_LUT_complete_cleaned.mat';

%% Load LUT
load(LUTname);

%% Plot
paramCount = size(LUT.table);
% For for each Dwake
for i = 1:paramCount(1)
    figure
    mesh(LUT.yWake, LUT.yaw, squeeze(LUT.table(i,2,:,:)))
    xlabel('yWake')
    ylabel('yaw')
    zlabel('DEL')
    title(['Dwake = ' num2str(LUT.Dwake(i)) ', Ufs = 8']);
    zlim([0 2E7]);
    xlim([-250 250]);
    ylim([-30 30]);
end
