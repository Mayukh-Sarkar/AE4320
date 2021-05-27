function [OLS] = OLS(C)
%{
the variable in the square brakets "ds" will be what the variable the
function outputs. The function is called using the file name "ds" and
inputting a variable that has at least 4 elements.
%}
%}
load('model.mat');
f = C(1).*omegaf.*exp(-omegaf*C(2))*C(3)^2./(omegaf.^2 + 2*C(4)*C(3).*omegaf + C(3)^2);
%{
The general form of a sine function. Be sure that all parentheses are
where they need to be. Alos, notice that the 1,2,3, and 4th constants are
referred to as C(1), C(2), C(3), and C(4).
%}

OLS = (f-Hp).^2;
%{
Defines an array where each element is equal to the difference squared
between corresponding elements of f and y. Notice the use of ".^2" instead
of "^2". The period tells Matlab to perform the calculation element by
element instead of through matrix multiplication.
%}

OLS = sum(OLS);
%{
Adds all elements of ds together. This last line will be the variable that
the function "ds" returns, and is the function to be minimized by
fminsearch.
%}