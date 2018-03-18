function best = jcbb(zz, H, best, level, ic, D_jc, Hx, S, Dz)

global State;
global Param;

	% if at leaf, return
	if level > size(zz,2)
		if nnz(H) > nnz(best)
			best = H;
		end

	else
		% prepare
		% at current level, check unassigned features, and create a "priority queue"(just use sorted array)
		nf = 0.5*(length(State.Ekf.mu) - 3); % total # of features
		q = [];
		f_idx = [];
		mu = State.Ekf.mu(1:3);
		% create queue
		for i = 1 : nf 
			if nnz(ismember(H, i))
				% already assigned
				continue;
			end

			% read individual compatibility from precomputed matrix
			q = [q, ic(level, i)];
			f_idx = [f_idx, i];
		end

		if isempty(q)
			% all features are assigned
			H(level:end) = 0;
			% best = jcbb(zz, H, best, level+1, ic, D_jc, Hx, S, Dz);
		else

			% descend sort
			[q, sorted_idx] = sort(q);
			f_idx = f_idx(sorted_idx); % sorted feature indices

			% do job
			allfail = false;
			% pruning if this is obviously worse than current best even if all of the rest measurements are able to be paired
			% might happens if there is a spurious at previous layers failing to pass every chi2test
			if nnz(H)+size(zz,2)-level+1 >= nnz(best) % still possible to be better than best

				for i = 1 : length(q)

					if q(i) > Param.ICthres
						% prune the rest
						allfail = true;
						break;
					end

					% incrementaly test JC
					[D_jc_tmp, Hx_tmp, S_tmp, Dz_tmp] = JC(zz(:,level), D_jc, Hx, S, Dz, f_idx(i));
					if D_jc_tmp < chi2inv(Param.JClevel, length(Dz_tmp))
						H(level) = f_idx(i);
						best = jcbb(zz, H, best, level+1, ic, D_jc_tmp, Hx_tmp, S_tmp, Dz_tmp);
						H(level) = 0;
					end
				end
			end

			if allfail || nnz(H)+size(zz,2)-level >= nnz(best) % must be sprious or possibly better if set spurious
				H(level) = 0;
				best = jcbb(zz, H, best, level+1, ic, D_jc, Hx, S, Dz);
			end
		end
	end
end