function q = facelift(p)
abs_p = abs(p);
k = sum(abs_p)/length(p)*10^-7;
for i = 1:length(p)
    if abs_p(i)<k
        q(i) = 0;
    else
        q(i) = p(i);
    end;
end;