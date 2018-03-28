# logistic3
Three Parameters logistic regression.</br>
One big holes into MatLab cftool function is the absence of Logistic Functions. In particular, The Three Parameters Logistic Regression or 3PL nonlinear regression model is commonly used for curve-fitting analysis in bioassays or immunoassays such as ELISA, RIA, IRMA or dose-response curves. It is characterized by it’s classic “S” or sigmoidal shape that fits the top plateaus of the curve, the EC50, and the slope factor (Hill's slope). This curve is symmetrical around its inflection point.
The 3PL equation is:
F(x) = D/(1+(x/C)^(-B))
where:<br/>
B = Hill's slope. The Hill's slope refers to the steepness of the curve. It could either be positive or negative.

C = Inflection point. The inflection point is defined as the point on the
curve where the curvature changes direction or signs. C is the concentration of analyte where y=D/2.

D = Maximum asymptote. In a bioassay where you have a standard curve, this can be thought of as the response value for infinite standard concentration.

In this submission there are 2 functions:
L3P - to find the 3 parameters and to fit your data (as calibrators...);
L3Pinv - to interpolate data of unknown samples onto calibrators curve.
Enjoy!

           Created by Giuseppe Cardillo
           giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2012) Three parameters logistic regression - There and back again
https://it.mathworks.com/matlabcentral/fileexchange/38124
