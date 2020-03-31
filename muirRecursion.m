clear all;
close all
syms z eq r;

n = 21;
rVal = 0.5;
numRecursions = 0;
[num, numRecursions] = A(n, z^(-1), r, numRecursions);
% den = A(n, z^(-1), -r);

coefficientsNumPre = fliplr(coeffs(expand((2/0.001)^0.5 * num * z^n), z));
% coefficientsDenPre = fliplr(coeffs(expand(den * z^n), z));

coefficientsNum = vpa(subs(coefficientsNumPre, r, rVal), 5)
% coefficientsDen = vpa(subs(coefficientsDenPre, r, rVal), 5)

function [eq, numRecursions] = A(n, z, r, numRecursions)
    numRecursions = numRecursions + 1;
    if n == 0
        eq = 1;
        return;
    end
    if mod(n, 2) == 0
        [eq, numRecursions] = A(n-1, 1/z, r, numRecursions);
    else
        [firstPart, numRecursions] = A(n-1, 1/z, r, numRecursions);
        [secondPart, numRecursions] = A(n-1, z, r, numRecursions);
        eq = firstPart - r / n * z^n * secondPart;
    end
end