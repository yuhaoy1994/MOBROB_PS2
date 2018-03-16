function Li = da_nn(z, R)
% perform nearest-neighbor data association

global Param;
global State;

Li = zeros(1,size(z,2));


% if no landmarks observed before
if length(State.Ekf.mu) <= 3
	for i = 1 : size(z,2)

		Li(i) = i;

	end
else

	mu = State.Ekf.mu(1:3);

	% from measurement to location
	lm = State.Ekf.mu(4:end);
	lm = reshape(lm, [2, length(lm)/2]);

	% measurement predicted from map
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
		
		delta_z = [z_pred(1,i) - z(1,:);
				  arrayfun(@wrapToPi, z_pred(2,i) - z(2,:))];
		costMat = [costMat; pdist2(zeros(1,2), delta_z', 'mahalanobis', S)];

	end


	while 1
		done = true;
		Li_tmp = [];
		[assign, cost] = munkres(costMat);
		% get assignment
		for i = 1 : size(assign,2)
			idx = find(assign(:,i));

			% check ambiguilty
			cost = costMat(:, i);
			if nnz(abs(cost - costMat(idx,i)) < 0.1) >= 2
				costMat(:,i) = []; % discard this measurement
				
				num_discard = nnz(Li == -1);
				Li(i+num_discard) = -1; % discarded measurement marked as -1
				done = false;
				break;
			end

			% check compartibility
			if costMat(idx,i)^2 < chi2inv(0.95,2)
				Li_tmp = [Li_tmp, idx];
			else
				% new landmark
				Li_tmp = [Li_tmp, size(lm,2)+1];
			end
		end		

		if done 
			break;
		end

	end

	% copy Li_tmp into Li
	j = 1;
	num_mapped = size(lm,2);
	for i = 1 : length(Li)
		if Li(i) == -1
			continue;
		else
			if Li_tmp(j) > size(lm,2)
				Li(i) = num_mapped + 1;
				num_mapped = num_mapped + 1;
			else
				Li(i) = Li_tmp(j);
			end
			j = j + 1;
		end
	end

end
