function x = L3Pinv(cf, y)
%L3PINV Inverse of the three-parameter logistic (3PL) equation.
%
%   x = L3Pinv(cf, y)
%
%   Description:
%       L3Pinv computes the inverse of the three-parameter logistic (3PL)
%       model used in L3P. Given a fitted 3PL curve or a numeric parameter
%       vector, and a set of response values y, this function returns the
%       corresponding x values that satisfy:
%
%           y = D / (1 + (x / C)^(-B))
%
%       where:
%           B : Hill slope (steepness and direction of the curve)
%           C : Inflection point (x at half-max, y = D/2)
%           D : Upper asymptote (maximum response)
%
%   Syntax:
%       x = L3Pinv(cf, y)
%
%   Inputs:
%       cf  - Either:
%              • A cfit object returned by L3P, containing the parameters
%                B, C, and D of the 3PL model, or
%              • A numeric vector [B, C, D] (1x3 or 3x1) with the 3PL
%                parameters.
%
%       y   - Numeric array of response values for which you want to
%             compute the corresponding x values. y must be real, finite,
%             and non-empty. It may be a scalar, vector, or matrix; the
%             output x will have the same size as y.
%
%   Outputs:
%       x   - Numeric array of x values such that the 3PL model evaluated
%             at x returns (approximately) the values in y. The size of x
%             matches the size of y.
%
%   Example:
%       % Fit a 3PL curve using L3P:
%       xdata = [0; 4.5; 10.6; 19.7; 40; 84; 210];
%       ydata = [0.0089; 0.0419; 0.0873; 0.2599; 0.7074; 1.528; 2.7739];
%
%       [cf, G] = L3P(xdata, ydata);
%
%       % Invert the fitted curve at a specific response:
%       yq = 1.782315;
%       xq = L3Pinv(cf, yq);
%
%       % Alternatively, using the explicit parameter vector:
%       params = coeffvalues(cf);   % [B, C, D]
%       xq2   = L3Pinv(params, yq); % should match xq
%
%   GitHub repository:
%       https://github.com/dnafinder/logistic3
%
%   Citation:
%       Cardillo G. (2025). logistic3: Three-parameter logistic regression
%       tools in MATLAB (L3P and L3Pinv). Available at:
%       https://github.com/dnafinder/logistic3
%
%   License:
%       This function is distributed under the terms specified in the
%       LICENSE file of the logistic3 repository.
%
%   Author:
%       Giuseppe Cardillo
%       giuseppe.cardillo.75@gmail.com
%
%   Created:
%       2012-01-01 (original concept)
%
%   Updated:
%       2025-11-19 (refactored and documented version)
%
%   Version:
%       1.1.0
%
%   See also:
%       L3P

% -----------------------------
% Input parsing and validation
% -----------------------------

p = inputParser;

% cf must be either a cfit-like object or a numeric vector of length 3.
addRequired(p, 'cf', @(v) ...
    isobject(v) || ...
    (isnumeric(v) && isvector(v) && numel(v) == 3 && ...
     all(isfinite(v(:)) & isreal(v(:)) & ~isnan(v(:)))));

% y must be numeric, real, finite, non-empty; any shape is allowed.
addRequired(p, 'y', @(v) ...
    validateattributes(v, {'numeric'}, ...
                       {'real', 'finite', 'nonnan', 'nonempty'}));

parse(p, cf, y);
cf = p.Results.cf;
y  = p.Results.y;

clear p;

% -----------------------------
% Extract 3PL parameters (B, C, D)
% -----------------------------
% If cf is a fit object, extract its coefficients; otherwise use the
% numeric vector directly.

if isobject(cf)
    params = coeffvalues(cf); % expects [B, C, D]
else
    % Ensure params is a row vector [B, C, D].
    params = cf(:).';
end

B = params(1);
C = params(2);
D = params(3);

% -----------------------------
% Basic sanity checks on y
% -----------------------------
% In the standard 3PL model, meaningful inversions usually require:
%   0 < y < D
% Values outside this range can lead to complex or non-physical x.
% We do not forbid these values, but we warn the user.
if any(y(:) <= 0)
    warning('L3Pinv:NonPositiveY', ...
        ['Some response values are <= 0. For the standard 3PL model, ' ...
         'meaningful inversion typically requires 0 < y < D.']);
end

if any(y(:) >= D)
    warning('L3Pinv:YSaturated', ...
        ['Some response values are >= D. For the standard 3PL model, ' ...
         'values at or above the upper asymptote may not invert cleanly.']);
end

% -----------------------------
% Invert the 3PL equation
% -----------------------------
% Starting from:
%   y = D / (1 + (x / C)^(-B))
%
% Rearranging:
%   D / y = 1 + (x / C)^(-B)
%   (x / C)^(-B) = D / y - 1
%   x / C       = (D / y - 1)^(-1 / B)
%   x           = C * (D / y - 1)^(-1 / B)
%
% The computation is performed element-wise over y.

ratio = (D ./ y) - 1;      % D/y - 1
x     = C .* (ratio .^ (-1 ./ B));

end
