function cdf = peak_thm_cdf(dimention, x, k, low, space)
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory 
% --dimention = 1 or 2
% --x is a vector
% --k is scalar
% --low is -infinity theoretically
% Optional 
% -- space is the spacing in cumtrapz for demention = 3
%--------------------------------------------------------------------------
% OUTPUT
% --cdf is cumulative distribution density, a vector with the same size as x
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
if dimention == 1
    den = @(x) peak_thm_density(1, x, k);
    cdf = 1:length(x);
    for i = 1:length(x)
        cdf(i) = integral(den, low, x(i));
    end
end

if dimention == 2
    den = @(x) peak_thm_density(2, x, k); 
    cdf = 1:length(x);
    for i = 1:length(x)
        cdf(i) = integral(den, low, x(i));
    end
end

if dimention == 3
    x_max = max(x);
    y = transpose(low:space:x_max);
    den = peak_thm_density(3, y, k); % try for dimention 1 and 2
    cdf_num = cumtrapz(y,den); % use trapz
    index = round((x-low)/space);% drop loop; check integral too
    cdf = cdf_num(index);
end


end