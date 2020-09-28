close all;
clear;
clc;

fname='inputParameters.txt';
nx=readParam(fname,'nx');
nz=readParam(fname,'nz');
Fs_x=readParam(fname,'Fs_x');
Fs_z=readParam(fname,'Fs_z');
gridSize=readParam(fname,'gridSize');
shotInterval=readParam(fname,'shotInterval');
travelTimeFileDir=readParam(fname,'travelTimeFileDir');
if exist(travelTimeFileDir, 'dir')
    rmdir(travelTimeFileDir,'s');
    mkdir(travelTimeFileDir);
else
    mkdir(travelTimeFileDir);
end
Fs=Fstar(Fs_z,Fs_x);
velocityModel=vmodel(nz,nx);

iz = Fs_z:nz+Fs_z-1;
ix = Fs_x:nx+Fs_x-1;
Fs_z2 = 2*Fs_z-1;
Fs_x2 = 2*Fs_x-1;
S = ones(nz+Fs_z2-1,nx+Fs_x2-1);
S(iz,ix) = 1./velocityModel;
S(nz+Fs_z,ix)=2*S(nz+Fs_z-1,ix)-S(nz+Fs_z-2,ix);
S(iz,nx+Fs_x)=2*S(iz,nx+Fs_x-1)-S(iz,nx+Fs_x-2);
S(nz+Fs_z,nx+Fs_x)=2*S(nz+Fs_z-1,nx+Fs_x-1)-S(nz+Fs_z-2,nx+Fs_x-2);

x = (0:nx) * gridSize;
z = (0:nz) * gridSize;

xShotGrid=1:shotInterval:nx;
nShots = length(xShotGrid);
xShot = (xShotGrid-1) * gridSize;

xRecGrid=1:nx;
xRec = (xRecGrid-1) * gridSize;

xShotAndRecGrid = union(xShotGrid, xRecGrid);
nShotsAndRecs = length(xShotAndRecGrid);

travelTime = zeros(nz, nx, nShotsAndRecs);

for ixs = 1:nShotsAndRecs 
    xs = xShotAndRecGrid(ixs); % shot position on x  
    travelTime(:, :, ixs) = Time2d(S,[1 xs], gridSize,nz,nx, Fs_z, Fs_x, Fs); % [1 xs]: shot position on [z, x]
    imagesc(x, z, travelTime(:, :, ixs));
    xlabel('Distance (m)'), ylabel('Depth (m)');
    title(sprintf('Travel time for Shot No. %d at x = %dm', ixs, x(xs)));
    drawnow
end

filenameTravelTime = [travelTimeFileDir, '/travelTime.mat'];
save(filenameTravelTime, 'travelTime', '-v7.3')


