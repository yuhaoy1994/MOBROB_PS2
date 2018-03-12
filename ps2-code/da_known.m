function Li = da_known(z)
% EKF-SLAM data association with known correspondences

global Param;
global State;

markerId = z;
Li = [];
for j = 1 : length(z)
	id = z(j);
	li = -1;
	for i = 1 : length(State.Ekf.da_known)

		if State.Ekf.da_known(i) == id;
			li = i;
			break;
		end

	end

	if li < 0
		li = length(State.Ekf.da_known) + 1;
		State.Ekf.da_known = [State.Ekf.da_known; markerId];
	end
	Li = [Li, li];
end