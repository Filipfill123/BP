function rotationAnimation(time_int,O_int,Quat_int,time_knot,O_knot,Quat_knot)

axesLength = max(max(O_knot))/8;

figure('Name','Lin seg. with poly. blend - animation')
ax = axes;

for i = 1:length(time_knot)
    R_knot{i} = Quaternion2rotMatrixAndAngularVel(Quat_knot(:,i));
end

for i = 1:length(time_int)
    R_int{i} = Quaternion2rotMatrixAndAngularVel(Quat_int(:,i));
end

for i = 1:length(time_knot)
    CSplot(ax,O_knot(:,i),R_knot{i},6,axesLength);
end
plot3(O_knot(1,:),O_knot(2,:),O_knot(3,:))
plot3(O_int(1,:),O_int(2,:),O_int(3,:),'Color','k','LineWidth',2)

xlabel('X')
xlabel('Y')
xlabel('Z')

axis equal

while(true)
p1 = [];
p2 = [];
p3 = [];
pause
    for i = 1:10:length(time_int)
        delete(p1)

        p1 = CSplot(ax,O_int(:,i),R_int{1,i},2,axesLength);
        pause(0.1)
    end
end

end
