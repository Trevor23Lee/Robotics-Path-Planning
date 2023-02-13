clear;
Arm = fopen("arm");
Arm_data = textscan(Arm, '%f %f %s');  %Reads Arm file
Offset_data = cell2mat(Arm_data(1));   % Gets the first column of the file, which is the Offest data
Angle_data = cell2mat(Arm_data(2));    % Gets the second column of the file, which is the Angle data
T = convertCharsToStrings(Arm_data(3)); % Gets the third column of the file, which is T(R or P)
T_data = T{1,1}; 

n = Offset_data(1,1);   % Reads and gets the n from the first column
lambda = Angle_data(1,1);  % Reads and gets the Gamma from the second column

Offset_data = Offset_data(2:end, 1);   % Gets the correct matrix(column) of Offset
Angle_data = Angle_data(2:end, 1);     % Gets the correct matrix(column) of Angle
T_data = T_data(2:end,1);

fclose(Arm);

Traj = dlmread("trajectory");    % Reads the trajectory file
m = Traj(1,1);                   % Gets the m from file
Traj = Traj(2:end, :);           % Gets the matrix of the trajectory file

sum_theta = 0;
x = 0;              % Variables for forward kinamatics
y = 0;

for i = 1 : n       % Loops through the number of joints
    sum_theta = Angle_data(i)+ sum_theta;   % Theta
    x = Offset_data(i) * cos(sum_theta) + x; % Gets the x-axis using the forward kinamatics equation
    y = Offset_data(i) * sin(sum_theta) + y; % Gets the y-axis using the forward kinamatics equation
end


[T_x,T_y] = size(T_data);
Rotations = "True";     % Initialize that it contains all rotations
for check = 1 : T_x     % Checks through out the T values to see if any of them are P(Prismatic Joints)
    if T_data(check,1) == "P"  % If so, then make Rotations false
        Rotations = "False";
    end
end

epsilon = 0.000001;     % Epsilon value variable
Identity = [1 0; 0 1];      % Identity matrix for Jacobian Lambda
Joint_ang = zeros(1,n);     % For final matrix
J = [0;0];      % Part of jacobian
theta = 0;
sum_theta = 0;    % Variables for Jacobian
length = 0;
Jacobian = zeros(2,1);  % For storing jacobian
max_delta = sum(Offset_data)*(0.01);  % Max dnorm delta 
boundary = sum(Offset_data);    % Boundary if all it is all rotations

clear_offset = Offset_data;
clear_angle = Angle_data;     % The clearing for Angle/Offset


for j = 1 : m
    actual_xy = [x,y];      % actual
    desired_xy = Traj(j,:);   % desired
    if Rotations == "True"    % If it did all contain rotations, then 
        if (norm(desired_xy)-boundary) > epsilon  % need to check the boundary
            p_closest_desired_x = desired_xy(1,1);  % if it is out then we need to 
            p_closest_desired_y = desired_xy(1,2);  % determine the closest reachable point
            txy = p_closest_desired_y/p_closest_desired_x;
            theta = atan(txy);                          % Give from TA, to help determine point
            desired_xy(1,1) = cos(theta) * sum(Offset_data);
            desired_xy(1,2) = sin(theta) * sum(Offset_data);
        end
    end
    error = ([desired_xy(1,1) - actual_xy(1,1), desired_xy(1,2) - actual_xy(1,2)]);

    while norm(error) > epsilon    % While error is greater than epsilon

         error = ([desired_xy(1,1) - actual_xy(1,1), desired_xy(1,2) - actual_xy(1,2)]);

         % Finds the correct step size - distance for a robot 
         % moving each iteration
         if norm(error) > max_delta
             Error_change = error/norm(error) * max_delta;
         else
             Error_change = error;
         end
         Jacobian = zeros(2,1);     % Zero out the Jacobian

         for joint = 1 : n
             J = zeros(2,1);
             if T_data(joint,1) == "R"          % If it is a rotary joint through the joints, getting Jacobian
                 for rot = joint : n
                     theta = Angle_data(1:rot,1);
                     sum_theta = sum(theta);
                     length = Offset_data(rot,1);
                     J(1,1) = length * (-sin(sum_theta)) + J(1,1);  % x value for rotary joints
                     J(2,1) = length * (cos(sum_theta)) + J(2,1);   % y value for rotary joints
                 end
             end


             if T_data(joint,1) == "P"          % If it is a prismatic joint through the joints, getting Jacobian
                 theta = Angle_data(1:joint,1);
                 sum_theta = sum(theta);
                 J(1,1) = cos(sum_theta);   % x value for prismatic joints
                 J(2,1) = sin(sum_theta);   % y value for prismatic joints
             end
             Jacobian = horzcat(Jacobian, J);   % Adding the segment J to the Jacobian
         end
         Jacobian = Jacobian(:,2:end);      % Gets rid of that first row of empty zeros
         
         Jacobian_Trans = Jacobian.';       % Gets the transpose of the Jacobian for equation below
         Jacobian_Lambda = Jacobian_Trans * inv((Jacobian*Jacobian_Trans) + (lambda^2 * Identity)); % Equation of J_lambda = J^T((J*J^T)+(lambda^2*I))^-1
         Error_change_Trans = transpose(Error_change);
         Angle_change = Jacobian_Lambda*Error_change_Trans;    % Computes the change in the angle

         % Checking if it will be a rotary of Prismatic joint then adds
         % into offset or angle
         for arp = 1 : n
             if T_data(arp) == "R"
                 Angle_data(arp,1) = Angle_data(arp,1) + Angle_change(arp,1);
             end
             if T_data(arp) == "P"
                 Offset_data(arp,1) = Offset_data(arp,1) + Angle_change(arp,1);
             end
         end

         % Calculating the new xy values using forward kinematics as
         % previously before
         x = 0;
         y = 0;
         sum_theta = 0;
         for d = 1 : n
             length = Offset_data(d,1);
             theta = Angle_data(1:d,1);
             sum_theta = sum(theta);
             x = x + length*cos(sum_theta);
             y = y + length*sin(sum_theta);
         end
         actual_xy = [x,y];
         error = ([desired_xy(1,1) - actual_xy(1,1), desired_xy(1,2) - actual_xy(1,2)]);
    end

    for cp = 1 : n                      % Checking if prismatic then angle is equal to offset
        if T_data(cp,1) == "P"
            Angle_data(cp,1) = Offset_data(cp,1);
        end
    end

    Joint_ang = vertcat(Joint_ang, transpose(Angle_data));  % Final matrix with angles/offsets
    Offset_data = clear_offset;     % Clearing offset and angle
    Angle_data = clear_angle;
end

Joint_ang = Joint_ang(2:end,:);   % Gets rid of the first row which is just zeros

[rows,cols] = size(Joint_ang);
file = fopen('angles','w');     % Open file 'angles' and be able to write into it
for k = 1 : rows                % Print out the Joint angle matrix which has the angles/offset
    fprintf(file, '\n');
    fprintf(file, '%f ', Joint_ang(k, :));
end