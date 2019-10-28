% interpolates (linearly) the values of x so that the nonzero entries of
% x span the range [a,b] 
function xx = linearize(x,a,b)
    xx = x;
    alpha = (b-a)./(max(max(x))-min(min(x(x>0))));
    beta = alpha*min(min(x(x>0)))-a;
    xx(x>0) = alpha.*x(x>0)-beta;
end

