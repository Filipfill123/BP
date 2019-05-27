function p = CSplot(ax,O,R,lineWidth,axisLength)
axes(ax)
hold on
p(1) = line([O(1);O(1)+axisLength*R(1,1)],[O(2);O(2)+axisLength*R(2,1)],[O(3);O(3)+axisLength*R(3,1)],'LineWidth',lineWidth,'Color','r');
p(2) = line([O(1);O(1)+axisLength*R(1,2)],[O(2);O(2)+axisLength*R(2,2)],[O(3);O(3)+axisLength*R(3,2)],'LineWidth',lineWidth,'Color','g');
p(3) = line([O(1);O(1)+axisLength*R(1,3)],[O(2);O(2)+axisLength*R(2,3)],[O(3);O(3)+axisLength*R(3,3)],'LineWidth',lineWidth,'Color','b');
