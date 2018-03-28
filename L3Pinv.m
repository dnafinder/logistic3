function x = L3Pinv(cf,y)
%L3PINV The inverse of the 3 parameters logistic equation.
% The Three Parameters Logistic Regression or 3PL nonlinear regression model
% is commonly used for curve-fitting analysis in bioassays or immunoassays
% such as ELISAs or dose-response curves. 
%
% Syntax: x=L3Pinv(cf,y)
% 
% Inputs: 
%           cf is the object containing the 3 parameters of logistic
%           equation computed by L3P function. Alternatively, it can be a
%           1x3 array.
%
%           y is the array of the response that you want to iterpolate. 
%
% Outputs:
%           x is the vector of interpolated data.
% 
% Example:
%
% xs=[0 4.5 10.6 19.7 40 84 210]; ys=[0.0089 0.0419 0.0873 0.2599 0.7074 1.528 2.7739];
%
% Calling on MatLab the function: [cf,G]=L3P(x,y);
%
% you will find the 5 parameters of this curve. 
% 
% Calling on MatLab the function L3Pinv(cf,1.782315);
% 
%           Answer is:
% ans =
%
%  100.0211
%
% Alternatively, you can do:
% 
% P=[1.5124 108.1762 3.7889]; L3Pinv(P,1.782315);
% 
% with the same result.
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% See also L3P, L4P, L4Pinv, L5P, L5Pinv
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2012) Three parameters logistic regression - There and back again
% 

%--------------------Input errors handling section-------------------------
p=inputParser;
addRequired(p,'cf',@(x) isobject(x) || (isvector(x) && length(x)==3 && all(isnumeric(x) & isreal(x) & isfinite(x) & ~isnan(x))))
addRequired(p,'y',@(x) validateattributes(x,{'numeric'},{'row','real','finite','nonnan','nonempty'}));
parse(p,cf,y)

if isobject(cf)
    p=coeffvalues(cf);
else
    p=cf;
end

%-------------------------Interpolate--------------------------------------
x=p(2).*((p(3)./y-1).^(-1/p(1)));
