function varargout = runsim(stepsOrData, pauseLen)

global Param;
global Data;
global State;
rng(5)
close all;
makemovie = false;


if ~exist('pauseLen','var')
    pauseLen = 0.3; % seconds
end

if makemovie
    try 
        votype = 'avifile';
        vo = avifule('video.avi', 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWritter';
        vo = VideoWritter('video', 'MPEG-4');
        set(vo, 'FrameRate', min(5, 1/pauseLen));
        open(vo); 
    end
end


% Initalize Params
%===================================================

Param.choice = 'sim';

% parameters for jcbb
Param.ICthres = chi2inv(0.99,2);
Param.JClevel = 0.99;

Param.initialStateMean = [180 50 0]';

% max number of landmark observations per timestep
Param.maxObs = 2;

% number of landmarks per sideline of field (minimum is 3)
Param.nLandmarksPerSide = 4;

% Motion noise (in odometry space, see p.134 in book).
Param.alphas = [0.05 0.001 0.05 0.01].^2; % std of noise proportional to alphas

% Standard deviation of Gaussian sensor noise (independent of distance)
Param.beta = [10, deg2rad(10)]; % [cm, rad]
Param.R = diag(Param.beta.^2);

% Step size between filter updates, can be less than 1.
Param.deltaT=0.1; % [s]

% threshold to determine whether an observation is new landmark or not
Param.nnthres = 130;


if isscalar(stepsOrData)
    % Generate a data set of motion and sensor info consistent with
    % noise models.
    numSteps = stepsOrData;
    Data = generateScript(Param.initialStateMean, numSteps, Param.maxObs, Param.alphas, Param.beta, Param.deltaT);
else
    % use a user supplied data set from a previous run
    Data = stepsOrData;
    numSteps = size(Data, 1);
    global FIELDINFO;
    FIELDINFO = getfieldinfo;
end
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = Param.initialStateMean;
State.Ekf.Sigma = zeros(3);
State.Ekf.da_known = [];

coeff_lx = []; coeff_cnt = 1;
coeff_rlx = [];
det_cov = [];
time_stamp = [];
gt_m = []; pred_m = []; cov_m = [];
for t = 1:numSteps
    plotsim(t);

    %=================================================
    % data available to your filter at this time step
    %=================================================
    u = getControl(t);
    z = getObservations(t);

    %=================================================
    %TODO: update your filter here based upon the
    %      motionCommand and observation
    %=================================================
    ekfpredict_sim(u);
    [a,b,c] = ekfupdate(z);
    gt_m = [gt_m, a]; pred_m = [pred_m, b]; cov_m = [cov_m, c]; % gather information for analyzing data association

    %=================================================
    %TODO: plot and evaluate filter results here
    %=================================================
    % plot robot 
    rob_pos = State.Ekf.mu(1:3);
    rob_pos_cov = State.Ekf.Sigma(1:3, 1:3);

    plotrobot( rob_pos(1), rob_pos(2), rob_pos(3), 'green', 1, 'green');
    plotcov2d( rob_pos(1), rob_pos(2), rob_pos_cov, 'blue', 0, 0, 0, 3);

    

    % plot landmark
    for i = 1 : (length(State.Ekf.mu)-3)/2
        mi = [State.Ekf.mu(3 + 2*i-1); State.Ekf.mu(3 + 2*i);];
        ci = State.Ekf.Sigma(3 + 2*i-1 : 3 + 2*i, 3 + 2*i-1 : 3 + 2*i); 
        plot(mi(1), mi(2), '*')
        plotcov2d( mi(1), mi(2), rob_pos_cov, 'blue', 0, 0, 0, 3);
    end


    drawnow;
    if pauseLen > 0
        pause(pauseLen);
    end

    % add frame
    if makemovie
        F = getframe(gcf);
        switch votype
        case 'avifile'
            vo = addframe(vo, F);
        case 'VideoWritter'
            writeVideo(vo, F);
        otherwise
            error('unrecognized type');
        end
    end


    % gather result for analyze covariance evolve
    if mod(t, 50) == 0
        l = 0.5*(length(State.Ekf.mu) - 3);
        if l >= 4
            C = State.Ekf.Sigma;
            coeff_lx= [coeff_lx, [C(4,6)/sqrt(C(4,4)*C(6,6)); C(4,8)/sqrt(C(4,4)*C(8,8)); C(6,8)/sqrt(C(6,6)*C(8,8))]]; coeff_cnt = coeff_cnt + 1;
            coeff_rlx = [coeff_rlx, C(1, [4,6,8])'./[sqrt(C(1,1)*C(4,4)); sqrt(C(1,1)*C(6,6)); sqrt(C(1,1)*C(8,8))]];
            det_cov = [det_cov, det(State.Ekf.Sigma)];
            time_stamp = [time_stamp, t];
        end
    end


end

if makemovie
    fprintf('Writing video ...');
    switch votype
    case 'avifile'
        vo = close(vo);
    case 'VideoWritter'
        close(vo);
    otherwise
        error('unrecognized type');
    end
    fprintf('done\n');
end


% plot task1 .4 
% figure;
% plot(time_stamp, coeff_lx(1,:),'r', time_stamp, coeff_lx(2,:), 'g', time_stamp, coeff_lx(3,:), 'b'); title('evolvement of correlation coefficient between 3 landmarks');
% saveas(gcf,'coeff_l.png');

% figure;
% plot(time_stamp, coeff_rlx(1,:),'r', time_stamp, coeff_rlx(2,:), 'g', time_stamp, coeff_rlx(3,:), 'b'); title('evolvement of correlation coefficient between robot and 3 landmarks');
% saveas(gcf,'coeff_lr.png');

% figure;
% plot(time_stamp, log(det_cov)); title('evolvement of determinate of covariance');
% saveas(gcf,'det.png');

% plot task2 .2
% figure;
% l = size(gt_m,2);
% scatter(1:l, gt_m(1,:), 'rx') , hold on; 
% scatter(1:l, pred_m(1,:), 'b');
% scatter(1:l, pred_m(1,:) + 3 * sqrt(cov_m(1,:)), 15, 'g', 'filled');
% scatter(1:l, pred_m(1,:) - 3 * sqrt(cov_m(1,:)), 15, 'g', 'filled'); legend('ground truth', 'nn estimate', '3 sigma bounds', '3 sigma bounds')
% title('DA NN: errors on x axis');
% saveas(gcf, 'da_nn_x.png');

% figure;
% l = size(gt_m,2);
% scatter(1:l, gt_m(2,:), 'rx') , hold on; 
% scatter(1:l, pred_m(2,:), 'b');
% scatter(1:l, pred_m(2,:) + 3 * sqrt(cov_m(2,:)), 15,'g' ,'filled');
% scatter(1:l, pred_m(2,:) - 3 * sqrt(cov_m(2,:)), 15, 'g', 'filled'); legend('ground truth', 'nn estimate', '3 sigma bounds', '3 sigma bounds')
% title('DA NN: errors on y axis');
% saveas(gcf, 'da_nn_y.png');


if nargout >= 1
    varargout{1} = Data;
end

%==========================================================================
function u = getControl(t)
global Data;
% noisefree control command
u = Data.noisefreeControl(:,t);  % 3x1 [drot1; dtrans; drot2]


%==========================================================================
function z = getObservations(t)
global Data;
% noisy observations
z = Data.realObservation(:,:,t); % 3xn [range; bearing; landmark id]
ii = find(~isnan(z(1,:)));
z = z(:,ii);

%==========================================================================
function plotsim(t)
global Data;

%--------------------------------------------------------------
% Graphics
%--------------------------------------------------------------

NOISEFREE_PATH_COL = 'green';
ACTUAL_PATH_COL = 'blue';

NOISEFREE_BEARING_COLOR = 'cyan';
OBSERVED_BEARING_COLOR = 'red';

GLOBAL_FIGURE = 1;

%=================================================
% data *not* available to your filter, i.e., known
% only by the simulator, useful for making error plots
%=================================================
% actual position (i.e., ground truth)
x = Data.Sim.realRobot(1,t);
y = Data.Sim.realRobot(2,t);
theta = Data.Sim.realRobot(3,t);

% real observation
observation = Data.realObservation(:,:,t);

% noisefree observation
noisefreeObservation = Data.Sim.noisefreeObservation(:,:,t);

%=================================================
% graphics
%=================================================
figure(GLOBAL_FIGURE); clf; hold on; plotfield(observation(3,:));

% draw actual path (i.e., ground truth)
plot(Data.Sim.realRobot(1,1:t), Data.Sim.realRobot(2,1:t), 'Color', ACTUAL_PATH_COL);
plotrobot( x, y, theta, 'black', 1, ACTUAL_PATH_COL);

% draw noise free motion command path
plot(Data.Sim.noisefreeRobot(1,1:t), Data.Sim.noisefreeRobot(2,1:t), 'Color', NOISEFREE_PATH_COL);
plot(Data.Sim.noisefreeRobot(1,t), Data.Sim.noisefreeRobot(2,t), '*', 'Color', NOISEFREE_PATH_COL);

for k=1:size(observation,2)
    rng = Data.Sim.noisefreeObservation(1,k,t);
    ang = Data.Sim.noisefreeObservation(2,k,t);
    noisy_rng = observation(1,k);
    noisy_ang = observation(2,k);

    % indicate observed range and angle relative to actual position
    plot([x x+cos(theta+noisy_ang)*noisy_rng], [y y+sin(theta+noisy_ang)*noisy_rng], 'Color', OBSERVED_BEARING_COLOR);

    % indicate ideal noise-free range and angle relative to actual position
    plot([x x+cos(theta+ang)*rng], [y y+sin(theta+ang)*rng], 'Color', NOISEFREE_BEARING_COLOR);
end
