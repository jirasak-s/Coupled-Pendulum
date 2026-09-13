function X = CoupledPendulum_Freq(wspan, wn, Phi, gamma)
    
    X = zeros(length(wspan), 2);

    for i = 1:length(wspan)
        E = zeros(2);
        for j = 1:2
            E(:,j) = Phi(:,j)*gamma(j)/(wn(j)^2 - wspan(i)^2);
        end
        X(i,:) = sum(E, 2);
    end
end