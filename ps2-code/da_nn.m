function Li = da_nn(z, R)
% perform nearest-neighbor data association

global Param;
global State;

Li = [];


% if no landmarks observed before
if length(State.Ekf.mu) <= 3
	for i = 1 : size(z,2)

		Li = [Li, i];

	end
else

	mu = State.Ekf.mu(1:3);

	% from measurement to location
	lm = State.Ekf.mu(4:end);
	lm = reshape(lm, [2, length(lm)/2]);

	% measurement prediction
	z_pred = [pdist2(mu(1:2)', lm');
			  arrayfun(@wrapToPi, atan2(lm(2,:)-mu(2), lm(1,:)-mu(1)) - mu(3))];

	% mahalanobis distance
	costMat = [];
	for i = 1 : size(lm,2)
		mi = State.Ekf.mu(3+2*i-1 : 3+2*i);
		% uncertainty estimation
		% Jacobian wrt rob pose & observed lm
		q = norm(mi - mu(1:2))^2;
		H = [sqrt(q)*(mu(1) - mi(1)), sqrt(q)*(mu(2) - mi(2)), 	 0, sqrt(q)*(mi(1)-mu(1)),     sqrt(q)*(mi(2)-mu(2));
			 (mi(2) - mu(2)), -(mi(1) - mu(1)), -q, -(mi(2) - mu(2)), (mi(1) - mu(1))] /q ; 

		A = zeros(5, length(State.Ekf.mu));
		A(1:3,1:3) = eye(3);
		A(4:5, 2*i+2:2*i+3) = eye(2);	
		% to state space
		H = H * A;
		% Compute Innovation Covariance
		S = H * State.Ekf.Sigma * H' + Param.R;	
		

		costMat = [costMat; pdist2(z_pred(:,i)', z', 'mahalanobis', S)];

	end

	[assign, cost] = munkres(costMat);

	% get assignment
	for i = 1 : size(z,2)
		idx = find(assign(:,i));
		% check if distance is within threshold
		if costMat(idx, i)^2 < 5.991 % 95% confidence
			Li = [Li, idx];	
		else
			% new landmark
			Li = [Li, size(lm,2)+1];
		end
	end

end
