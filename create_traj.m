function[time_length, cesta] = create_traj(trajectory)

time_length = length(trajectory.time);
cesta = [];
for i = 1:time_length
    cesta = [cesta,trajectory.signals.values(:,:,i)];
end

end