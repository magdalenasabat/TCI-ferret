function newRGB = myAdjustBrightness(rgb, scalar)
    % Ensure input is a valid RGB triplet
    if any(rgb < 0) || any(rgb > 1) || length(rgb) ~= 3
        error('RGB must be a 1x3 vector with values between 0 and 1');
    end

    if scalar < 0
        error('Scalar must be non-negative');
    end

    if scalar < 1
        % Darken: scale down the RGB values
        newRGB = rgb * scalar;
    else
        % Lighten: interpolate towards white
        newRGB = rgb + (1 - rgb) * (scalar - 1);
    end

    % Clamp values to [0, 1]
    newRGB = max(min(newRGB, 1), 0);
end
