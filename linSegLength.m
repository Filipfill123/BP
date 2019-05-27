function arcLength = linSegLength(Pvect,Ovect,Ivect);
% lin segments length
% pouze z prvnich 3 souradnic (x,y,z)
numOfPoints = size(Pvect,2);

for i = 1:numOfPoints - 1
    if i == 1
        arcLength(i) = norm(Ivect(1:3,i)-Pvect(1:3,i));
    elseif i == numOfPoints - 1
        arcLength(i) = norm(Ovect(1:3,i-1)-Pvect(1:3,i+1));
    else
        arcLength(i) = norm(Ivect(1:3,i)-Ovect(1:3,i-1));
    end
end
