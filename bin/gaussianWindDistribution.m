function [ Udir_weight_vector ] = gaussianWindDistribution( windDirectionsRange , plotResults)
    x = [-8+min(windDirectionsRange):0.1:max(windDirectionsRange)+8];

    Udir_weight_vector = zeros(length(windDirectionsRange),1);
    for Udir_vector_step = 1:length(windDirectionsRange)
        Udir_weight_vector(Udir_vector_step,1) = normpdf(windDirectionsRange(Udir_vector_step),0,8);
    end
    Udir_weight_sum = sum(Udir_weight_vector);

    beta = (1/Udir_weight_sum);
    norm = beta*normpdf(x,0,8);
    Udir_weight_vector = zeros(length(windDirectionsRange),1);
    for Udir_vector_step = 1:length(windDirectionsRange)
        Udir_weight_vector(Udir_vector_step,1) = beta*normpdf(windDirectionsRange(Udir_vector_step),0,8);
    end
    Udir_weight_sum = sum(Udir_weight_vector);
    
    % Figure
    if plotResults
        grid on
        hold on
        plot(x,norm)
        plot((windDirectionsRange),Udir_weight_vector,'r.','MarkerSize',20)
        title('Normal distribution of wind directions')
        xlabel('Degrees [^o]')
        ylabel('Probability [-]')
    end
end

