function [arcLength] = polySegLength(u_knot,polySegPar);
% poly segments length
for i = 1:length(polySegPar)
    fun = @(u)integrandPoly(u,polySegPar{i});
    j = 2*i;
    arcLength(i) = quad(fun,u_knot(j),u_knot(j+1));
end


