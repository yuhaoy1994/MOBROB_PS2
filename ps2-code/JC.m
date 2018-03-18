function [D_jc, Hx, S, Dz] = JC(zz, D_jc, Hx, S, Dz, fidx)

global State;
global Param;
	
	mu = State.Ekf.mu(1:3);
	mi = [State.Ekf.mu(3+2*fidx-1); State.Ekf.mu(3+2*fidx)];
	% Jacobian
	q = norm(mi - mu(1:2))^2;
	Hi = [sqrt(q)*(mu(1) - mi(1)), sqrt(q)*(mu(2) - mi(2)), 0 , sqrt(q)*(mi(1)-mu(1)), sqrt(q)*(mi(2)-mu(2));
		 (mi(2) - mu(2))		 , -(mi(1) - mu(1))		  , -q, -(mi(2) - mu(2))	 , (mi(1) - mu(1))      ] /q ; 

	A = zeros(5, length(State.Ekf.mu));
	A(1:3,1:3) = eye(3);
	A(4:5, 2*fidx+2:2*fidx+3) = eye(2);	
	% to state space
	Hi = Hi * A;
	% innovation
	Sij = Hi * State.Ekf.Sigma * Hi' + Param.R;	

	% meassurement prediction
	switch Param.choice
	case 'sim'
		z_pred = [norm(mi-mu(1:2));
				  wrapToPi(atan2(mi(2)-mu(2), mi(1)-mu(1)) - mu(3))];
	case 'vp'
		z_pred = [norm(mi-mu(1:2));
				  wrapToPi(atan2(mi(2)-mu(2), mi(1)-mu(1)) - mu(3) + pi/2)];
	end

	delta_z = [z_pred(1) - zz(1);
			   wrapToPi(z_pred(2) - zz(2))];

	if isempty(Hx)
		Hx = Hi;
		S = Sij;
		Dz = delta_z;
		D_jc = delta_z' / S * delta_z;
	else
		% increntally evaluate JC following Jose Neira & Juan D. Tardos, 2001
		wi  = Hi(:,1:3) * State.Ekf.Sigma(1:3,1:3) * Hx(:,1:3)' + Hi(:,4:end) * State.Ekf.Sigma(4:end,4:end) * Hx(:,4:end)';
		inv_Ni = Sij  - wi/S*wi';
		Li = -inv_Ni\wi/S;
		try
			D_jc = D_jc + Dz'*Li'*inv_Ni*Li*Dz + 2*delta_z'*Li*Dz + delta_z'/inv_Ni*delta_z;
		catch ME
			keyboard;
			rethrow(ME);
		end
		S = [S, wi';wi, Sij];
		Hx = [Hx; Hi];
		Dz = [Dz; delta_z];
	end

end