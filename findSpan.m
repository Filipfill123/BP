function seg = findSpan(t,t_vect)
% nalezeni segmentu v intervalu t_vect, kde se nachazi hodnota t

[minVal,index] = min(abs(t_vect - t));

if t < t_vect(index)
    seg = index - 1;
else
    if t >= t_vect(end)
        seg = index-1;
    else
        seg = index;
    end
end
