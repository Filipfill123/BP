function T = DH(DHpar)

d = DHpar(1);
theta = DHpar(2);
a = DHpar(3);
alpha = DHpar(4);

Trans_z = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, d; 0, 0, 0, 1];
Rot_z = [cos(theta), -sin(theta), 0, 0; sin(theta), cos(theta), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
Trans_x = [1, 0, 0, a; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
Rot_x = [1, 0, 0, 0; 0, cos(alpha), -sin(alpha), 0; 0, sin(alpha), cos(alpha), 0;  0, 0, 0, 1];

T = Trans_z*Rot_z*Trans_x*Rot_x;

