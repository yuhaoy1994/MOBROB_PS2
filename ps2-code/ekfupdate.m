function ekfupdate(z)
% EKF-SLAM update step for both simulator and Victoria Park data set

global Param;
global State;


% do pairing 
% returns state vector indices pairing observations with landmarks
switch lower(Param.dataAssociation)
    case 'known'
        Li = da_known(z(3,:));
    case 'nn'
        Li = da_nn(z(1:2,:), Param.R);
    case 'jcbb'
        Li = da_jcbb(z(1:2,:), Param.R);
    otherwise
        error('unrecognized data association method: "%s"', Param.dataAssociation);
end

for i = 1 : size(z,2)
	MarkerId = z(3,i);



	% see if we already known this marker


	if 2*Li(i) > length(State.Ekf.mu) - 3
		initialize_new_landmark(z(:,i), Param.R);
	end
	% get old estimation of ith landmark
	mi = [State.Ekf.mu(3+2*Li(i)-1); State.Ekf.mu(3+2*Li(i))];
	% current prediction of rob pose
	mu = State.Ekf.mu(1:3);

	% predict measurement
	z_hat = [norm(mi - mu(1:2)); wrapToPi(atan2(-mu(2)+mi(2), -mu(1)+mi(1)) - mu(3))];

	% Jacobian wrt rob pose & observed lm
	q = norm(mi - mu(1:2))^2;
	H = [sqrt(q)*(mu(1) - mi(1)), sqrt(q)*(mu(2) - mi(2)), 	 0, sqrt(q)*(mi(1)-mu(1)),     sqrt(q)*(mi(2)-mu(2));
		 (mi(2) - mu(2)), -(mi(1) - mu(1)), -q, -(mi(2) - mu(2)), (mi(1) - mu(1))] /q ; 

	A = zeros(5, length(State.Ekf.mu));
	A(1:3,1:3) = eye(3);
	A(4:5, 2*Li(i)+2:2*Li(i)+3) = eye(2);	
	% to state space
	H = H * A;

	% Compute Innovation Covariance
	S = H * State.Ekf.Sigma * H' + Param.R;

	% Kalman gain
	K = State.Ekf.Sigma * H' / S;

	% Correction
	State.Ekf.mu = State.Ekf.mu + K*[z(1,i) - z_hat(1); wrapToPi(z(2,i) - z_hat(2))];
	State.Ekf.Sigma = (eye(length(State.Ekf.mu)) - K*H) * State.Ekf.Sigma;



end


