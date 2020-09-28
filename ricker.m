function [rw,t] = ricker(f,n,dt,t0,t1)

switch nargin
    case 0
        f  = 20;
        n  = 100;
        dt = 0.001;
        t0 = 1/f;
        is2d = false;
    case 1
        n = 100;
        dt = 0.001;
        t0 = 1/f;
        is2d = false;
    case 2
        dt = 0.001;
        t0 = 1/f;
        is2d = false;
    case 3
        t0 = 1/f;
        is2d = false;
    case 4 % use all values
        is2d = false;
    case 5 % use all inputs
        is2d = true;
    otherwise
        warning('RICKER:tooManyInputs','Ignoring additional inputs')
end

% Create the wavelet and shift in time if needed
T = dt*(n-1);
t = 0:dt:T;
tau = t-t0;
if ~is2d
    s = (1-2*tau.*tau*f^2*pi^2).*exp(-tau.^2*pi^2*f^2);
else
    [t1,t2] = meshgrid(tau,t-t1);
    s = (1-2*(t1.^2+t2.^2)*f^2*pi^2).*exp(-(t1.^2+t2.^2)*pi^2*f^2);
end

if nargout == 0
        plot(t,s)
        xlabel('Time (s)')
        ylabel('Amplitude')
        title('Ricker Wavelet')
else
        rw = s;
end