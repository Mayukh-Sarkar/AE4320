function [fitresult, gof] = createFit2(Freq_fit, sef_fit)
%CREATEFIT2(FREQ_FIT,SEF_FIT)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : Freq_fit
%      Y Output: sef_fit
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 24-May-2021 11:03:36


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( Freq_fit, sef_fit );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.101900586811651;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
set(gca, 'YScale', 'log')
legend( h, 'sef_fit vs. Freq_fit', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Freq_fit', 'Interpreter', 'none' );
ylabel( 'sef_fit', 'Interpreter', 'none' );
grid on

