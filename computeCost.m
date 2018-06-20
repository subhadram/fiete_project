function J = computeCost(X, y, theta)
%COMPUTECOST Compute cost for linear regression
%   J = COMPUTECOST(X, y, theta) computes the cost of using theta as the
%   parameter for linear regression to fit the data points in X and y

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 


% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost of a particular choice of theta
%               You should set J to the cost.
sm = 0.0;

for n = 1:m
    sm = sm + (theta(1) + theta(2)*X(n,2) -y(n))^2.0; 
    %((theta(1) + theta(2)*X(n)) -
end

J = sm/(2*m);



% =========================================================================

end
