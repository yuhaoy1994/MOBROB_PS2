function initialize_new_landmark(z, R)

global Param;
global State;

mu = State.Ekf.mu;
r = z(1);
theta = z(2) + mu(3);

% predict average location
mi_x = mu(1) + r * cos(theta);
mi_y = mu(2) + r * sin(theta);

% compute uncertainty 
% Jacobian wrt rob state
Fx = [1, 0, -r*sin(theta);
	  0, 1,  r*cos(theta)];

% Jacobian wrt senser noise
Fe = [cos(theta), -sin(theta);
	  sin(theta), cos(theta)];

% uncertainty
C = Fx * State.Ekf.Sigma(1:3,1:3) * Fx' + Fe * R * Fe';

% Augment state
Sigma = zeros(length(State.Ekf.mu)+2, length(State.Ekf.mu)+2);
Sigma(1:length(State.Ekf.mu), 1:length(State.Ekf.mu)) = State.Ekf.Sigma;
Sigma(length(State.Ekf.mu)+1:end, length(State.Ekf.mu)+1:end) = C;
State.Ekf.Sigma = Sigma; 

State.Ekf.mu = [State.Ekf.mu; mi_x; mi_y];