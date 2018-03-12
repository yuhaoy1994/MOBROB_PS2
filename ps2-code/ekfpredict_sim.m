function ekfpredict_sim(u)
% EKF-SLAM prediction for simulator process model

global Param;
global State;

% compute odometry noise
R = diag([Param.alphas(1)*u(1)^2 + Param.alphas(2)*u(2)^2;
		  Param.alphas(3)*u(2)^2 + Param.alphas(4)*(u(1)^2 + u(3)^2);
		  Param.alphas(1)*u(3)^2 + Param.alphas(2)*u(2)^2]);

% predict rob state mean
rob_mu = State.Ekf.mu(1:3);
th = rob_mu(3) + u(1);
State.Ekf.mu(1:3) = rob_mu + [u(2) * cos(th); u(2) * sin(th); u(1)+u(3)];

% compute jacobian
% Jacobian wrt state
Gx = eye(3);
Gx(1,3) = -u(2)*sin(th);
Gx(2,3) = u(2)*cos(th);

% Jacobian wrt noise
Ge = [-u(2)*sin(th), cos(th), 0;
 	   u(2)*cos(th), sin(th), 0;
 	   1		   ,	   0, 1];


% predict covariance
b11 = Gx * State.Ekf.Sigma(1:3,1:3) * Gx' + Ge * R * Ge';
b12 = Gx * State.Ekf.Sigma(1:3,4:end); % uncorrelated noise assumption
b21 = b12';

State.Ekf.Sigma(1:3,1:3) = b11;
State.Ekf.Sigma(1:3,4:end) = b12;
State.Ekf.Sigma(4:end,1:3) = b21;
