function best = da_jcbb(z, R)
% perform joint-compatability branch and bound data association

global Param;
global State;

if length(State.Ekf.mu) < 4
	best = zeros(1, size(z,2));
else
	mu = State.Ekf.mu(1:3);

	% get mapped landmarks
	lm = State.Ekf.mu(4:end);
	lm = reshape(lm, [2, length(lm)/2]);

	% measurement predicted from map
	switch Param.choice
	case 'sim'
		z_pred = [pdist2(mu(1:2)', lm');
				  arrayfun(@wrapToPi, atan2(lm(2,:)-mu(2), lm(1,:)-mu(1)) - mu(3))];
	case 'vp'
		z_pred = [pdist2(mu(1:2)', lm');
				  arrayfun(@wrapToPi, atan2(lm(2,:)-mu(2), lm(1,:)-mu(1)) - mu(3) + pi/2)];
	end

	% compute mahalanobis distance between each pair of meassurements and features
	ic = [];
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
		
		try 
			ic = [ic, pdist2(delta_z', zeros(1,2), 'mahalanobis', S)];
			% costMat = [costMat; diag(delta_z'/S*delta_z)'];
		catch ME
	        if strcmp(ME.identifier, 'stats:pdist2:InvalidCov') % sometimes I get this wield error from MATLAB when S is PSD
	            [~, p] = chol(S);
	            if p == 0 && rank(S) == size(S,1)
	                ic = [ic, diag(delta_z'/S*delta_z)];
	            else
	                rethrow(ME);
	            end
	        else
	            rethrow(ME);
	        end
		end
	end

	level = 1;
	H = zeros(1, size(z,2)); % hypothesis; if spurious fill 0
	best = [];
	D_jc = 0;
	Hx = [];
	S = [];
	Dz = [];
	best = jcbb(z, H, best, level, ic, D_jc, Hx, S, Dz);
end

% init all spurious as newlandmark
n_newlm = 1;
for i = 1 : length(best)
	if best(i) == 0
		best(i) = 0.5*(length(State.Ekf.mu) - 3) + n_newlm;
		n_newlm = n_newlm + 1;
	end
end

