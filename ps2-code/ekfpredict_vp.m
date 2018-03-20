function ekfpredict_vp(u, dt)
% EKF-SLAM prediction for Victoria Park process model

global Param;
global State;


ve = u(1);
al = u(2);
th = State.Ekf.mu(3);

a = Param.a; b = Param.b; L = Param.L; H = Param.H;

vc = ve / (1 - Param.H/Param.L * tan(al));

% predict mean state
dx = dt * [vc*cos(th) - vc/L*tan(al)*(a*sin(th) + b*cos(th));
		   vc*sin(th) + vc/L*tan(al)*(a*cos(th) - b*sin(th));
		   vc*tan(al)/L];

State.Ekf.mu(1:3) = State.Ekf.mu(1:3) + dx;

% predict covariance
% jacobian wrt state
Gx = eye(3);
Gx(1,3) = -dt * (vc*sin(th) + vc*tan(al)*(a*cos(th) - b*sin(th))/L);
Gx(2,3) = dt * (vc*cos(th) - vc*tan(al)*(a*sin(th) + b*cos(th))/L);
% jacobian wrt control noise
Ge = dt * [cos(th)-tan(al)/L*(a*sin(th)+b*cos(th)), -vc/L/(cos(al)^2)*(a*sin(th)+b*cos(th));
		   sin(th)+tan(al)/L*(a*cos(th)-b*sin(th)), vc/L/(cos(al)^2)*(a*cos(th)-b*sin(th));
		   tan(al)/L 							  , vc/L/(cos(al)^2)];
b11 = Gx * State.Ekf.Sigma(1:3,1:3) * Gx' + Ge * Param.Qu * Ge' + Param.Qf;
b12 = Gx * State.Ekf.Sigma(1:3,4:end);
b21 = b12';

State.Ekf.Sigma(1:3,1:3) = b11;
State.Ekf.Sigma(1:3,4:end) = b12;
State.Ekf.Sigma(4:end,1:3) = b21;
