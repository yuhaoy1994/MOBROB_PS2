function runvp(nSteps,pauseLen)
close all
global Param;
global State;
global Data;

if ~exist('nSteps','var') || isempty(nSteps)
    nSteps = inf;
end

if ~exist('pauseLen','var')
    pauseLen = 0; % seconds
end

Data = load_vp_si();

% Initalize Params
%===================================================

Param.choice = 'vp';
Param.ICthres = chi2inv(0.99,2);
Param.JClevel = 0.99;
% vehicle geometry
Param.a = 3.78; % [m]
Param.b = 0.50; % [m]
Param.L = 2.83; % [m]
Param.H = 0.76; % [m]

% 2x2 process noise on control input
sigma.vc = 0.02; % [m/s]
sigma.alpha = 2*pi/180; % [rad]
Param.Qu = diag([sigma.vc, sigma.alpha].^2);

% 3x3 process noise on model error
sigma.x = 0.1; % [m]
sigma.y = 0.1; % [m]
sigma.phi = 0.5*pi/180; % [rad]
Param.Qf = diag([sigma.x, sigma.y, sigma.phi].^2);

% 2x2 observation noise
sigma.r = 0.05; % [m]
sigma.beta = 1*pi/180; % [rad]
Param.R = diag([sigma.r, sigma.beta].^2);
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = [Data.Gps.x(2), Data.Gps.y(2), 36*pi/180]';
State.Ekf.Sigma = zeros(3);


global AAr;
AAr = [0:360]*pi/360;


% figure(1); clf;
% axis equal;


traj = [];
ci = 1; % control index
t = min(Data.Laser.time(1), Data.Control.time(1));
for k=1:min(nSteps, length(Data.Laser.time))
    
    while (Data.Control.time(ci) < Data.Laser.time(k))
       % control available
       dt = Data.Control.time(ci) - t;
       t = Data.Control.time(ci);
       u = [Data.Control.ve(ci), Data.Control.alpha(ci)]';

       ekfpredict_vp_sol(u, dt);
       traj = [traj, [State.Ekf.mu(1); State.Ekf.mu(2)]];
       
       ci = ci+1;
    end
    
    % observation available
    dt = Data.Laser.time(k) - t;
    t = Data.Laser.time(k);
    z = detectTreesI16(Data.Laser.ranges(k,:));

    

    ekfupdate(z);
    traj = [traj, [State.Ekf.mu(1); State.Ekf.mu(2)]];

    for i = 1 : length(Data.Gps.time)
      if Data.Gps.time(i) > t
        break;
      end
    end

    if i >= length(Data.Gps.time)
      warning('Gps used up\n');
    end

    gt_traj = [Data.Gps.x(1:i)'; Data.Gps.y(1:i)'];

    doGraphics(z, traj, gt_traj);
    drawnow;
    if pauseLen > 0
        pause(pauseLen);
    end
end




%==========================================================================
function doGraphics(z,traj, gt_traj)
% Put whatever graphics you want here for visualization
%
% WARNING: this slows down your process time, so use sparingly when trying
% to crunch the whole data set!

global Param;
global State;

% plot the robot and 3-sigma covariance ellipsoid
plotbot(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.mu(3), 'black', 1, 'blue', 1);
hold on;

plotcov2d( State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma(1:2,1:2), 'blue', 0, 'blue', 0, 3);


% restrict view to a bounding box around the current pose
% BB=20;
% axis([[-BB,BB]+State.Ekf.mu(1), [-BB,BB]+State.Ekf.mu(2)]);
axis auto

% project raw sensor detections in global frame using estimate pose
xr = State.Ekf.mu(1);
yr = State.Ekf.mu(2);
tr = State.Ekf.mu(3);
for k=1:size(z,2)
    r = z(1,k);
    b = z(2,k);
    xl = xr + r*cos(b+tr-pi/2);
    yl = yr + r*sin(b+tr-pi/2);
    plot([xr; xl], [yr; yl],'m',xl,yl,'m*');
end


for i = 1 : (length(State.Ekf.mu) - 3)/2

  plotcov2d( State.Ekf.mu(3+2*i-1), State.Ekf.mu(3+2*i), State.Ekf.Sigma(3+2*i-1:3+2*i, 3+2*i-1:3+2*i), 'red', 0, 'red', 0, 3 );

end

% plot trajectory
plot(traj(1,:), traj(2,:), 'r'), hold on
plot(gt_traj(1,:), gt_traj(2,:), 'b');

hold off;

