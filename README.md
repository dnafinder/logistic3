📘 Overview
logistic3 is a lightweight MATLAB toolbox that implements the three-parameter logistic (3PL) regression model, widely used in bioassays and immunoassays such as ELISA, RIA, IRMA, and dose-response curves. The 3PL model captures the sigmoidal (“S-shaped”) relationship between analyte concentration and measured response, including the upper asymptote (D), slope factor (B), and inflection point (C, EC50).

The mathematical form of the 3PL equation is:
    F(x) = D / (1 + (x / C)^(-B))

where:
- B is the Hill slope (curve steepness and direction; positive or negative)
- C is the inflection point (x value at half-max response, y = D/2)
- D is the maximum asymptote (upper plateau of the curve)

This repository provides two functions:
- L3P    : fits the 3PL model to your data
- L3Pinv : evaluates the inverse of the 3PL model (e.g., estimating concentrations from measured responses)

✨ Features
- Implements canonical 3PL curve fitting for biological and analytical assays
- Accepts single-measurement data or replicate matrices (with automatic weighting)
- Uses nonlinear least squares regression with adjustable starting values and parameter bounds
- Provides goodness-of-fit statistics: SSE, R², adjusted R², degrees of freedom, RMSE
- Supports forward evaluation (model prediction) and inverse evaluation (value interpolation)
- Based entirely on MATLAB built-in tools plus Curve Fitting Toolbox

📥 Installation
1. Download the logistic3 repository from GitHub:
   https://github.com/dnafinder/logistic3

2. Place the folder in any directory of your choice.

3. Add the folder to your MATLAB path:
      addpath('path_to_logistic3')

4. Verify MATLAB can find the functions:
      which L3P
      which L3Pinv

⚙️ Requirements
- MATLAB (recent releases recommended)
- Curve Fitting Toolbox (for fit, fittype, fitoptions, and cfit objects)

📈 Usage
Fitting a 3PL model with L3P:

    x = [0; 4.5; 10.6; 19.7; 40; 84; 210];
    y = [0.0089; 0.0419; 0.0873; 0.2599; 0.7074; 1.528; 2.7739];

    [cf, G] = L3P(x, y);

Plotting the result:

    plot(x, y, 'ro');
    hold on;
    plot(cf, 'r');
    hold off;

Interpolating unknown samples with L3Pinv:

    response = 1.8;
    x_est = L3Pinv(cf, response);

Alternatively, using an explicit parameter vector:

    params = [B, C, D];  % your 3PL parameters
    x_est = L3Pinv(params, response);

🔢 Inputs
L3P:
- x  : Column vector of concentrations (N×1)
- y  : Column vector of responses OR matrix of replicates (N×M)
- st : Optional starting values [B0 C0 D0]
- L  : Optional lower bounds
- U  : Optional upper bounds

L3Pinv:
- cf : cfit object returned by L3P OR numeric vector [B C D]
- y  : Query response values, any shape (scalar, vector, matrix)

📤 Outputs
L3P returns:
- cf : cfit object representing the fitted 3PL model
- G  : Structure of goodness-of-fit metrics (SSE, R², adjR², RMSE)

L3Pinv returns:
- x : Numeric array of interpolated x values with the same shape as y

🧠 Interpretation
- B > 0 : increasing sigmoidal curve
- B < 0 : decreasing sigmoidal curve
- C     : EC50 or midpoint of the curve
- D     : upper asymptote (plateau)

A good fit typically shows:
- High R² and adjusted R²
- Low SSE and RMSE
- A fitted curve visually matching the data points

📌 Notes
- When providing a replicate matrix, L3P automatically averages responses and uses standard deviations as weights.
- Ensure that the measured responses y fall within a meaningful range for 3PL inversion (usually 0 < y < D).
- Nonlinear regression can be sensitive to poor starting values or unrealistic bounds.

🧾 Citation
If you use logistic3 in publications or analyses, please cite:

Cardillo G. (2025). logistic3: Three-parameter logistic regression tools in MATLAB (L3P and L3Pinv). Available at:
https://github.com/dnafinder/logistic3

👤 Author
Giuseppe Cardillo  
Email: giuseppe.cardillo.75@gmail.com  
GitHub: https://github.com/dnafinder

📄 License
This project is distributed under the terms specified in the LICENSE file available at:
https://github.com/dnafinder/logistic3
