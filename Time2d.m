function T = Time2d(S,Shot,dx, nz,nx,Fs_z, Fs_x, Fs)
%2D time field
T0 = 1.e8;
Fs_z2 = 2*Fs_z-1;
Fs_x2 = 2*Fs_x-1;
zs = Shot(1) + Fs_z - 1;
xs = Shot(2) + Fs_x - 1;


mxV = max(S(:));

T = ones(nz+Fs_z2,nx+Fs_x2) * T0;
M = T;

%%
M(Fs_z:nz+Fs_z,Fs_x:nx+Fs_x) = 0;
iz = Fs_z:nz+Fs_z-1;
ix = Fs_x:nx+Fs_x-1;

T(zs,xs) = 0;
M(zs,xs) = T0;

z1 = -Fs_z+1:Fs_z-1; z2 = -Fs_z+1:Fs_z-2; z3 = z1+zs;
x1 = -Fs_x+1:Fs_x-1; x2 = -Fs_x+1:Fs_x-2; x3 = x1+xs;
AS = S(z2+zs,x2+xs);
TT = T(z3,x3);
T(z3,x3) = min(reshape(Fs*AS(:)+T(zs,xs),Fs_z2,Fs_x2),TT);
mxT = max(max(T(zs-1:zs+1,xs-1:xs+1)));

while true
    indx = T+M <= mxT+mxV;
    
    if isempty(indx)
         indx = M == 0;
    end
    [idz,idx] = find(indx);
    M(indx) = T0;
    for i = 1:length(idz)
        z = idz(i);
        x = idx(i);
        mxT = max(mxT,T(z,x));
        AS = S(z+z2,x+x2);
        z3 = z+z1;
        x3 = x+x1;
        TT = T(z3,x3);
        T(z3,x3) = min(reshape(Fs*AS(:)+T(z,x),Fs_z2,Fs_x2),TT);
    end %for
    if prod(double(all(M(iz,ix))))
        break;
    end
    mxT = max(max(T(idz,idx)));
end %while

T = T(iz,ix)*dx;