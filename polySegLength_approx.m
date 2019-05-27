function [arcLength] = polySegLength_approx(Pvect,Ovect,Ivect);
% poly segments length Approx.(normalized segments u = <0,1>)
for i = 1:size(Ovect,2);
    arcLength(i) = norm(Ivect(:,i)-Pvect(:,i+1)) + norm(Ovect(:,i)-Pvect(:,i+1));
end
