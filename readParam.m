function out=readParam(fname,paramName)
fid=fopen(fname);
cac = textscan( fid,'%s%f%s');
var=cac{1,1};
value=cac{1,2};
vs=cac{1,3};
n=length(var);
n1=length(value);
for k=1:n
    if(strcmp(paramName,var(k)))
        out=value(k);
        if(out==0)
            out=char(vs(k));
        end
        break;
    end
end

