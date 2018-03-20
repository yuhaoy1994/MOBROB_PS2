function [gt_m, pred_m, cov_m] = ekfupdate(z)
% EKF-SLAM update step for both simulator and Victoria Park data set

global Param;
global State;
global FIELDINFO;

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

switch Param.updateMethod
case 'batch'
	m = [];
	mu = State.Ekf.mu(1:3);

	for i = 1 : size(z,2)
		if 2*Li(i) > length(State.Ekf.mu) - 3
			initialize_new_landmark(z(:,i), Param.R);
		end
	end

	H = [];
	for i = 1 : size(z,2)

		if Li(i) == -1
			continue;
		end

		% get landmark location
		mi = [State.Ekf.mu(3+2*Li(i)-1); State.Ekf.mu(3+2*Li(i))];
		m = [m, mi];
		% Jacobian
		q = norm(mi - mu(1:2))^2;
		Hi = [sqrt(q)*(mu(1) - mi(1)), sqrt(q)*(mu(2) - mi(2)), 	 0, sqrt(q)*(mi(1)-mu(1)), sqrt(q)*(mi(2)-mu(2));
			 (mi(2) - mu(2)), -(mi(1) - mu(1)), -q, -(mi(2) - mu(2)), (mi(1) - mu(1))] /q ; 

		A = zeros(5, length(State.Ekf.mu));
		A(1:3,1:3) = eye(3);
		A(4:5, 2*Li(i)+2:2*Li(i)+3) = eye(2);	
		% to state space
		H = [H; Hi * A];		 
	end
	S = H * State.Ekf.Sigma * H' + kron(eye(nnz(Li ~= -1)), Param.R);
	K = State.Ekf.Sigma * H' / S;

	% current prediction of state
	mu = repmat(mu, [1, size(m,2)]);

	% predict measurement
	switch Param.choice
	case 'sim'
		z_hat = [ sqrt(sum((m - mu(1:2,:)).^2));
			  arrayfun(@wrapToPi, atan2(m(2,:) - mu(2,:), m(1,:) - mu(1,:)) - mu(3,:));
			];
	case 'vp'
		z_hat = [ sqrt(sum((m - mu(1:2,:)).^2));
			  arrayfun(@wrapToPi, atan2(m(2,:) - mu(2,:), m(1,:) - mu(1,:)) - mu(3,:) + pi/2);
			];
	end



	dz = z(1:2,Li~=-1) - z_hat;
	dz(2,:) = arrayfun(@wrapToPi, dz(2,:));

	% reshape to column vector
	dz = reshape(dz, [nnz(Li~=-1)*2, 1]);    
    % correction
	State.Ekf.mu = State.Ekf.mu + K * dz;
	State.Ekf.Sigma = (eye(length(State.Ekf.mu)) - K*H) * State.Ekf.Sigma;



case 'seq'
	gt_m = [];
	pred_m = [];
	cov_m = [];
	for i = 1 : size(z,2)
		try 
		if Li(i) == -1
			continue;
		end
	catch
	keyboard 
	end

		% see if we already known this marker
		if 2*Li(i) > length(State.Ekf.mu) - 3
			initialize_new_landmark(z(:,i), Param.R);
        end

        if strcmp(Param.choice, 'sim')
        	gt_m = [gt_m, [FIELDINFO.MARKER_X_POS(z(3,i)); FIELDINFO.MARKER_Y_POS(z(3,i))]];
        	pred_m = [pred_m, [State.Ekf.mu(3+2*Li(i)-1); State.Ekf.mu(3+2*Li(i))]];
        	cov_m = [cov_m, [State.Ekf.Sigma(3+2*Li(i)-1,3+2*Li(i)-1); State.Ekf.Sigma(3+2*Li(i),3+2*Li(i))]];
    	end
        try 
		% get old estimation of ith landmark
		mi = [State.Ekf.mu(3+2*Li(i)-1); State.Ekf.mu(3+2*Li(i))];
        catch ME
            save('err.mat');
            rethrow(ME);
        end
		% current prediction of rob pose
		mu = State.Ekf.mu(1:3);

		% predict measurement
		switch Param.choice
		case 'sim'
			z_hat = [norm(mi - mu(1:2)); wrapToPi(atan2(-mu(2)+mi(2), -mu(1)+mi(1)) - mu(3))];
		case 'vp'
			z_hat = [norm(mi - mu(1:2)); wrapToPi(atan2(-mu(2)+mi(2), -mu(1)+mi(1)) - mu(3) + pi/2)];
		end

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

end





