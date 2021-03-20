function [S, SHalf, SBar, addPointsAt] = setTube(N, NnonExtended, n, setToOnes)
    if setToOnes
        addPointsAt = ceil(0.5 * (N+1));
        S = ones(N+1, 1);
    else
        lengths = [0.708, 0.177, 0.711, 0.306, 0.254, 0.502];
        radii = [0.0069, 0.0072, 0.0069, 0.0071, 0.0075, 0.0107]; % two radii for tuning slide

        lengthN = round(NnonExtended * lengths ./ sum(lengths));
        addPointsAt = round(lengthN(1) + lengthN(2) * 0.5) + (N-NnonExtended) * 0.5; % indicate split of two connected schemes (including offset if N differs from NnonExtended

        inner1 = ones(lengthN(1), 1) * radii(1);
        inner2 = ones(lengthN(3), 1) * radii(3);
        gooseneck = ones(lengthN(4), 1) * radii(4);
        tuning = linspace(radii(5), radii(6), lengthN(5))';

        x0 = 0.0174; 
        b = 0.0063;
        flare = 0.7;
        bellL = lengthN(end);

        bell = b * ((lengths(6):-lengths(6) / (bellL - 1):0) + x0).^(-flare);


    %     pointsLeft = N - length([mp, m2t, bell]);
    %     tube = linspace(m2t(end), m2t(end), pointsLeft);    % tube
        totLengthN = length(inner1) + length(inner2) + length(gooseneck) + length(tuning) + length(bell);
    %     lengthN(2) = lengthN(2) + (N - NnonExtended + 1);
        slide = ones(N - totLengthN, 1) * radii(2);
        totRadii = [inner1; slide; inner2; gooseneck; tuning; bell'];

        % True geometry
        S = totRadii.^2 * pi;
    end
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;                  % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                        % mu_{x-}S_{l+1/2}
end