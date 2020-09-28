close all;
clear;
clc;
EPSILON = 1e-3;


fname='inputParameters.txt';
nx=readParam(fname,'nx');
nz=readParam(fname,'nz');
gridSize=readParam(fname,'gridSize');
dx=gridSize;dz=gridSize;
% Define frequency parameter for ricker wavelet
f = readParam(fname,'f');
shotWatchInterval=readParam(fname,'shotWatchInterval');
shotInterval=readParam(fname,'shotInterval');
shotGatherFileDir=readParam(fname,'shotGatherFileDir');
if ~exist(shotGatherFileDir, 'dir')
    error('Step1_modeling: Data needs prepared!');
end
travelTimeFileDir=readParam(fname,'travelTimeFileDir');
if ~exist(travelTimeFileDir, 'dir')
    error('Step2_time2d: Traveltime needs prepared!');
end
velocityModel=vmodel(nz,nx);
nBoundary = 20;
x = (0:nx) * dx;
z = (0:nz) * dz;
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));
dt = 0.75*(dz/vmax/sqrt(2));
nt = round(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1);
t  = (0:nt-1).*dt;
% grids and positions of shot array
xShotGrid=1:shotInterval:nx;
nShots = length(xShotGrid);
xShot = (xShotGrid-1) * dx;
xRecGrid=1:nx;
xRec = (xRecGrid-1) * gridSize;

xShotAndRecGrid = union(xShotGrid, xRecGrid);
nShotsAndRecs = length(xShotAndRecGrid);

filenameTravelTime = [travelTimeFileDir, '/travelTime.mat'];
load(filenameTravelTime);
Stacked = zeros(nz, nx);
%% Process Shots - Kirchhoff Migration
Ssignal=ricker(f, nt, dt);
[~,kstart]=max(Ssignal);
hShotPos = plot(xShot(1), z(1), 'w*');
subplot(2,2,1);
imagesc(x, z, velocityModel)
xlabel('Distance (m)'); ylabel('Depth (m)');
title('Velocity Model');
hold on;
hShotPos = plot(xShot(1), z(1), 'w*');
hold off;
colormap(seismic);
for ixs = 1:nShots
    xs = xShotGrid(ixs); % shot position on x
    shotGatherFilename = [shotGatherFileDir, sprintf('/shotGather%d.mat', ixs)];
    load(shotGatherFilename); 
     %data = data(nBoundary+1:end-nBoundary, :);
    tic;
     M = migrate(travelTime, xRecGrid, xShotAndRecGrid, data(kstart:end,:), dt, nz, xs);
    timeKirchoffMigration = toc;
    fprintf('Kirchoff Migration for Shot No. %d at x = %d, elapsed time = %fs\n', xs, x(xs), timeKirchoffMigration);
    set(hShotPos, 'XData', x(xs), 'YData', z(1)); 
    Stacked = Stacked + M;
    
    subplot(2,2,2)
    imagesc(x, z, Stacked)
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Stacked Image');
    %colormap(seismic);
    % caxis([-135 135])
    
    subplot(2,2,3)
    imagesc(x,t,data)
    xlabel('Distance (m)'); ylabel('Time (s)');
    title(sprintf('Current Shot No. %d at x = %dm', ixs, x(xs)));
    caxis([-0.1 0.1])
    
    subplot(2,2,4)
    imagesc(x,t,M)
    xlabel('Distance (m)'); ylabel('Time (s)');
    title(sprintf('Current Migrated Shot No. %d at x = %dm', ixs, x(xs)));
    %colormap(seismic);  
    %set(hShotPos, 'XData', x(xs));
    
    drawnow

end
