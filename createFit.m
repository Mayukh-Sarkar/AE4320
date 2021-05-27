function [fitresult, gof] = createFit(omegaf, magh)
%CREATEFIT(OMEGAF,MAGH)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : omegaf
%      Y Output: magh
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 25-May-2021 13:10:53


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( omegaf, magh );

% Set up fittype and options.
ft = fittype( 'p1.*x*exp(-x*p2)*p3^2./(x^2 + 2*p4*p3.*x + p3^2);', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.141886338627215 0.421761282626275 0.915735525189067 0.792207329559554];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'magh vs. omegaf', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'omegaf', 'Interpreter', 'none' );
ylabel( 'magh', 'Interpreter', 'none' );
grid on

