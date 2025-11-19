function [cf, G] = L3P(x, y, varargin)
%L3P Three-parameter logistic regression (3PL).
%
%   [cf, G] = L3P(x, y)
%   [cf, G] = L3P(x, y, st)
%   [cf, G] = L3P(x, y, st, L, U)
%
%   Description:
%       L3P fits a three-parameter logistic (3PL) model to experimental
%       data, typically arising from dose-response or standard curves in
%       bioassays (e.g., ELISA). The 3PL model captures the sigmoidal shape
%       of the response, including the upper asymptote, the inflection
%       point (EC50), and the slope (Hill slope). The 3PL model used is:
%
%           F(x) = D / (1 + (x / C)^(-B))
%
%       where:
%           B : Hill slope (steepness and direction of the curve).
%           C : Inflection point (x value where y = D/2).
%           D : Upper asymptote (maximum response).
%
%   Inputs:
%       x   - Column vector (N x 1) of x-values (e.g., concentrations).
%       y   - Either:
%               • Column vector (N x 1) of responses, or
%               • Matrix (N x M) of replicate responses for each x.
%             In the matrix case, L3P computes the mean response per row,
%             and uses the row-wise standard deviation as weights.
%
%       st  - (Optional) 1x3 row vector of initial parameter estimates:
%                 [B0, C0, D0]
%             If omitted or empty, reasonable starting values are estimated
%             from (x, y).
%
%       L   - (Optional) 1x3 row vector of lower bounds for [B, C, D].
%             If omitted or empty, bounds are inferred from the data:
%                 • Default L = [0, 0, 0], except:
%                 • If the initial slope is negative, L(1) = -Inf.
%
%       U   - (Optional) 1x3 row vector of upper bounds for [B, C, D].
%             If omitted or empty, bounds are inferred from the data:
%                 • Default U = [Inf, Inf, Inf], except:
%                 • If the initial slope is negative, U(1) = 0.
%
%   Outputs:
%       cf  - A cfit object representing the fitted 3PL model:
%                 cf(x) = D / (1 + (x / C)^(-B))
%
%       G   - Structure with goodness-of-fit statistics, including:
%                 G.sse        : sum of squares due to error
%                 G.rsquare    : coefficient of determination (R^2)
%                 G.dfe        : degrees of freedom in the error
%                 G.adjrsquare : adjusted R^2
%                 G.rmse       : root mean square error
%
%   Example:
%       x = [0; 4.5; 10.6; 19.7; 40; 84; 210];
%       y = [0.0089; 0.0419; 0.0873; 0.2599; 0.7074; 1.528; 2.7739];
%
%       [cf, G] = L3P(x, y);
%
%       % Plot data and fitted curve
%       figure;
%       hold on;
%       plot(x, y, 'ro', 'DisplayName', 'Data');
%       plot(cf, 'r');
%       legend('show', 'Location', 'best');
%       hold off;
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
%       L3Pinv

% -----------------------------
% Input parsing and validation
% -----------------------------

% Use inputParser to validate mandatory and optional inputs.
p = inputParser;

% x must be a numeric, real, finite, non-empty column vector.
addRequired(p, 'x', @(v) validateattributes(v, ...
    {'numeric'}, {'column', 'real', 'finite', 'nonnan', 'nonempty'}));

% y must be a numeric, real, finite, non-empty 2D array (vector or matrix).
addRequired(p, 'y', @(v) validateattributes(v, ...
    {'numeric'}, {'2d', 'real', 'finite', 'nonnan', 'nonempty'}));

% Optional starting values for [B, C, D].
addOptional(p, 'st', [], @(v) isempty(v) || ...
    validateattributes(v, {'numeric'}, {'row', 'real', 'finite', 'nonnan', 'numel', 3}));

% Optional lower bounds for [B, C, D].
addOptional(p, 'L', [], @(v) isempty(v) || ...
    validateattributes(v, {'numeric'}, {'row', 'real', 'nonnan', 'numel', 3}));

% Optional upper bounds for [B, C, D].
addOptional(p, 'U', [], @(v) isempty(v) || ...
    validateattributes(v, {'numeric'}, {'row', 'real', 'nonnan', 'numel', 3}));

% Parse all inputs
parse(p, x, y, varargin{:});

x  = p.Results.x;
y  = p.Results.y;
st = p.Results.st;
L  = p.Results.L;
U  = p.Results.U;

clear p;

% Ensure that x and y are compatible in size.
assert(size(x, 1) == size(y, 1), ...
    'L3P:SizeMismatch', 'x and y must have the same number of rows.');
assert(size(x, 1) >= 3, ...
    'L3P:NotEnoughData', 'Not enough data points. At least 3 rows are required.');

% -----------------------------
% Preprocessing of y (replicates)
% -----------------------------

% If y is a matrix, compute row-wise means and standard deviations.
if ~isvector(y)
    % Standard deviation for each row (used as weights).
    we = std(y, 0, 2);
    % Mean response for each x.
    y  = mean(y, 2);
else
    % Ensure y is a column vector and set weights to zero (no weighting).
    y  = y(:);
    we = zeros(size(x));
end

% -----------------------------
% Initial parameter estimates
% -----------------------------
% B (Hill slope) is approximated with the slope of the line joining the
% first and last data points.
% C (inflection point) is approximated as the x at which y is closest to
% half of the maximum observed response.
% D (upper asymptote) is approximated as max(y).

% Basic slope estimate from first to last point.
slope = (y(end) - y(1)) / (x(end) - x(1));

% If required, infer starting values from the data.
if isempty(st)
    % Approximate inflection point as the concentration where response
    % is closest to half of the maximum response.
    yMid     = max(y) / 2;
    [~, Idx] = min(abs(y - yMid));
    
    % Use the sign of the slope for B, x(Idx) for C, and max(y) for D.
    st = [sign(slope), x(Idx), max(y)];
end

% -----------------------------
% Parameter bounds
% -----------------------------
% Default logic:
%   - C and D are constrained to be non-negative by default.
%   - If the slope is positive, B is constrained to be non-negative.
%   - If the slope is negative, B is allowed to be negative, but not
%     positive (upper bound = 0).

if isempty(L)
    L = zeros(1, 3); % [B_lower, C_lower, D_lower]
    if slope < 0
        % Allow negative slopes if data suggest a decreasing curve.
        L(1) = -Inf;
    end
end

if isempty(U)
    U = Inf(1, 3); % [B_upper, C_upper, D_upper]
    if slope < 0
        % If the curve is decreasing, prevent positive B.
        U(1) = 0;
    end
end

% -----------------------------
% Define fit options and model
% -----------------------------
% Configure nonlinear least squares fitting with starting points and bounds.
fo = fitoptions('Method', 'NonlinearLeastSquares', ...
                'StartPoint', st, ...
                'Lower', L, ...
                'Upper', U);

% If weights are available (i.e., y was a matrix and all std > 0),
% include them in the fitting process.
if all(we)
    set(fo, 'Weights', we);
end

% Define the three-parameter logistic model as a fittype object.
ft = fittype('D / (1 + (x / C)^(-B))', ...
    'dependent',  {'y'}, ...
    'independent',{'x'}, ...
    'coefficients',{'B', 'C', 'D'});

% -----------------------------
% Perform the fit
% -----------------------------
% Use the Curve Fitting Toolbox function FIT to estimate the parameters.
[cf, G] = fit(x, y, ft, fo);

end
