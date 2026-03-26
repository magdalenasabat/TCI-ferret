function v = robust_variance(X, f)

% subtract median
X = X - repmat(median(X), size(X,1), 1);

% estimate standard deviation from central samples
sd = (quantile(X, 0.5+f/2)-quantile(X, 0.5-f/2)) ...
    / (norminv(0.5+f/2,0,1) - norminv(0.5-f/2,0,1));

v = sd.^2;