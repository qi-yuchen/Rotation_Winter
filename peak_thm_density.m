function den = peak_thm_density(dimention, x, k)
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory 
% --dimention = 1, 2 or 3
% --x is a vector
% --k is scalar
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% --den is distribution density, a vector with the same size as x
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
if dimention == 1
    den = sqrt(3-k^2)/sqrt(6*pi)*exp(-3*x.^2/(2*(3-k^2)))...
    + 2*k*x*sqrt(pi)/sqrt(6).*normpdf(x)...
    .*normcdf(k*x/sqrt(3-k^2));
end

if dimention == 2
    den = sqrt(3)*k^2*(x.^2-1).*normpdf(x)...
    .*normcdf(k*x/sqrt(2-k^2)) + k*x*sqrt(3*(2-k^2))/(2*pi).*...
    exp(-x.^2/(2-k^2)) + sqrt(6)/sqrt(pi*(3-k^2))...
    *exp(-3*x.^2/(2*(3-k^2)))...
    .*normcdf(k*x/sqrt(3-k^2)/sqrt(2-k^2));
end

if dimention == 3
    a = k^2*((1-k^2)^3 + 6*(1-k^2)^2 + 12*(1-k^2) + 24)...
        /(4*(3-k^2)^2);
    b = (2*(1-k^2)^3 + 3*(1-k^2)^2 + 6*(1-k^2))/(4*(3-k^2));
    c = 7-k^2+(1-k^2)*(3*(1-k^2)^2+12*(1-k^2)+28)/2/(3-k^2);
    sigma1 = [3/2,-1;-1,(3-k^2)/2];
    s = mvncdf([zeros(length(x), 1),k*x/sqrt(2)],[],sigma1);
    sigma2 = [3/2,-1/2;-1/2,(2-k^2)/2];
    t = mvncdf([zeros(length(x), 1),k*x/sqrt(2)],[],sigma2);
    lat_part = (a*x.^2 + b + 3/2).*exp(-k^2*x.^2/(2*(3-k^2)))...
        /sqrt(2*(3-k^2)).*normcdf(2*k*x/sqrt(3-k^2)/sqrt(5-3*k^2))...
        + (k^2*(2-k^2)/4*x.^2 - k^2*(1-k^2)/2 - 1).*...
        exp(-k^2*x.^2/(2*(2-k^2)))/sqrt(2*(2-k^2)).*...
        normcdf(k*x/sqrt(2-k^2)/sqrt(5-3*k^2))...
        + c*k*x.*exp(-3*k^2*x.^2/(2*(5-3*k^2)))/4/sqrt(pi)...
        /(3-k^2)/sqrt(5-3*k^2)...
        + sqrt(pi)*k^3/4*x.*(x.^2-3).*(s + t);
    den = 144*normpdf(x)/(29*sqrt(6)-36).*lat_part;
end


end
