function initialize_new_landmark(z, R)

global Param;
global State;

mu = State.Ekf.mu;
l = length(mu);
r = z(1);
theta = z(2) + mu(3);

% predict average location
switch Param.choice
	case 'sim'
		mi_x = mu(1) + r * cos(theta);
		mi_y = mu(2) + r * sin(theta);
		% Jacobian wrt rob state
		Fx = [1, 0, -r*sin(theta);
			  0, 1,  r*cos(theta)];

		% Jacobian wrt senser noise
		Fe = [cos(theta), -r*sin(theta);
			  sin(theta), r*cos(theta)];
	case 'vp'
		mi_x = mu(1) + r * cos(theta - pi/2);
		mi_y = mu(2) + r * sin(theta - pi/2);
		% Jacobian wrt rob state
		Fx = [1, 0, -r*sin(theta - pi/2);
			  0, 1,  r*cos(theta - pi/2)];

		% Jacobian wrt senser noise
		Fe = [cos(theta-pi/2), -r*sin(theta-pi/2);
			  sin(theta-pi/2), r*cos(theta-pi/2)];
end
		

% compute uncertainty 


% augment covariance
cxx = State.Ekf.Sigma(1:3,1:3);

cxi = Fx * cxx;
cii = Fx * cxx * Fx' + Fe * R * Fe';


Sigma = zeros(l+2, l+2);
Sigma(1:l, 1:l) = State.Ekf.Sigma;
Sigma(l+1:end, 1:3) = cxi;
Sigma(1:3, l+1:end) = cxi';

if l > 3
	cxj = State.Ekf.Sigma(1:3, 4:end);
	cij = Fx * cxj;
	Sigma(l+1:end, 4:l) = cij;
	Sigma(4:l, l+1:end) = cij';
end

Sigma(l+1:end, l+1:end) = cii;
State.Ekf.Sigma = Sigma;

State.Ekf.mu = [State.Ekf.mu; mi_x; mi_y];