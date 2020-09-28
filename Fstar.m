function A=Fstar(sz,sx)
sz2 = 2*sz-1;
sx2 = 2*sx-1;
sz1=sz2-1;
sx1=sx2-1;
nrow=sz2*sx2;
ncol=sz1*sx1;
A=zeros(nrow,ncol);
nray=0;
rayxz=zeros(2,1000);
temp=zeros(2,1000);

for kz=1:sz2
    z0=kz-1;
    for kx=1:sx2
        x0=kx-1;
        nray=nray+1;
        dxx=(sx-kx);
        dzz=(sz-kz);
        if(dxx==0 || dzz==0)
            if(dxx==0 && dzz~=0)
                if(dzz>0)
                    np=0;
                    for kk=kz:sz
                        np=np+1;
                        temp(1,np)=x0;
                        temp(2,np)=kk-1;
                    end
                else
                    np=0;
                    for kk=kz:-1:sz
                        np=np+1;
                        temp(1,np)=x0;
                        temp(2,np)=kk-1;
                    end  
                end
             else
                if(dxx>0)
                    np=0;
                    for kk=kx:sx
                        np=np+1;
                        temp(1,np)=kk-1;
                        temp(2,np)=z0;
                    end                  
                else
                    np=0;
                    for kk=kx:-1:sx
                        np=np+1;
                        temp(1,np)=kk-1;
                        temp(2,np)=z0;
                    end          
                end
            end           
        else
            slop=dzz/dxx;
            rslop=1./slop;
            if(slop>0)
                if(dxx>0)      % I
                    seg=1;
                    rayxz(1,seg)=x0;
                    rayxz(2,seg)=z0;
                    x=x0;
                    for ix=kx+1:sx
                        seg=seg+1;
                        x=x+1;
                        z=slop*(x-x0);
                        rayxz(1,seg)=x;
                        rayxz(2,seg)=z+z0;  
                    end
                    z=z0; 
                    for iz=kz+1:sz
                        seg=seg+1;
                        z=z+1;
                        x=(z-z0)*rslop;
                        rayxz(1,seg)=x+x0;
                        rayxz(2,seg)=z;  
                    end
                else            % IV
                    seg=1;
                    rayxz(1,seg)=x0;
                    rayxz(2,seg)=z0;
                    x=x0;
                    for ix=kx-1:-1:sx
                        seg=seg+1;
                        x=x-1;
                        z=slop*(x-x0);
                        rayxz(1,seg)=x;
                        rayxz(2,seg)=z+z0;  
                    end
                    z=z0; 
                    for iz=kz-1:-1:sz
                        seg=seg+1;
                        z=z-1;
                        x=(z-z0)*rslop;
                        rayxz(1,seg)=x+x0;
                        rayxz(2,seg)=z;  
                    end
                end
            else
                if(dxx<0)    %II
                    seg=1;
                    rayxz(1,seg)=x0;
                    rayxz(2,seg)=z0;
                    x=x0;
                    for ix=kx-1:-1:sx
                        seg=seg+1;
                        x=x-1;
                        z=slop*(x-x0);
                        rayxz(1,seg)=x;
                        rayxz(2,seg)=z+z0;  
                    end
                    z=z0; 
                    for iz=kz+1:sz
                        seg=seg+1;
                        z=z+1;
                        x=(z-z0)*rslop;
                        rayxz(1,seg)=x+x0;
                        rayxz(2,seg)=z;  
                    end
                else                 %III
                    seg=1;
                    rayxz(1,seg)=x0;
                    rayxz(2,seg)=z0;
                    x=x0;
                    for ix=kx+1:sx
                        seg=seg+1;
                        x=x+1;
                        z=slop*(x-x0);
                        rayxz(1,seg)=x;
                        rayxz(2,seg)=z+z0;  
                    end
                    z=z0; 
                    for iz=kz-1:-1:sz
                        seg=seg+1;
                        z=z-1;
                        x=(z-z0)*rslop;
                        rayxz(1,seg)=x+x0;
                        rayxz(2,seg)=z;  
                    end
                end
            end
            %sorting
            [out,idx] = sort(rayxz(1,1:seg));
            rayxz(1,1:seg)=out;
            rayxz(2,1:seg)=rayxz(2,idx);
            temp(:,1)=rayxz(:,1);
            np=1;
            for k=2:seg
                aa=rayxz(:,k)-rayxz(:,k-1);
                dist=sqrt(aa(1)*aa(1)+aa(2)*aa(2));
                if(dist>1.e-5)
                    np=np+1;
                    temp(:,np)=rayxz(:,k);
                end
            end
        end
        for k=2:np
            aa=temp(:,k)-temp(:,k-1);
            dist=sqrt(aa(1)*aa(1)+aa(2)*aa(2));
            aa=0.5*(temp(:,k)+temp(:,k-1));
            indx=floor(aa(1));
            indz=floor(aa(2));
            ind=indz*sx1+indx+1;
            A(nray,ind)=dist;
        end      
    end
end
            
        
