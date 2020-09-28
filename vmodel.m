function v=vmodel(nz,nx)
v=zeros(nz,nx);
nz1=floor(nz/5);
v(1:nz1,:)=2000;
v(nz1+1:end,:)=3500;