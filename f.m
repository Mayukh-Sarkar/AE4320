function y = f(C,omegaf)
y = C(1).*omegaf.*exp(-omegaf*C(2))*C(3)^2./(omegaf.^2 + 2*C(4)*C(3).*omegaf+ C(3)^2);

end