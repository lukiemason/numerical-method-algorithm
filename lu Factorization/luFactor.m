function[L,U,P] = luFactor(A)
%This function uses lu factorization along with partial pivoting to
%determine the upper and lower triangle matrix, when given a square matrix.
%   
%Inputs:
% A - the coefficient matrix
%Outputs:
% L - Lower triangle
% U - Upper Triangle
% P - Pivot Matrix
%
% Author: Luke Mason
% Date: 3/26/2018
%
[n,m] = size(A); %This allows us to determine the size of the matrix
if m~=n %makes sure the matrix is a square
    error('input matrix must be a square')
end
if nargin > 1 %right number of inputs
    error('too many inputs, silly(;')
end
[n,m] = size(A);
L = eye(n); %identity matrix
P = eye(n); %identity matrix
U = A
for i = 1: m-1
    Pivot = max(abs(U(i:m,i))); %record how the pivoting takes place
    for j = i:m
        if (abs(U(j,i)) == Pivot)
            independent = j;
        end
    end
    U([i,independent],i:m) = U([independent,i],i:m); %changing rows Upper
    L([i,independent],1:i-1) = L([independent,i],1:i-1); %changing rows Lower
    P([i,independent],:) = P([independent,i],:); %changing rows pivot
    for j = i+1:m
        L(j,i) = U(j,i)/U(i,i); %sovling for lower triangle
        U(j,i:m) = U(j,i:m) - L(j,i)*U(i,i:m); %solvong for upper triangle
    end
end
display(L) %
display(U) %displaying answeres
display(P) %
end

