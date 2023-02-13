% Plot the desired end effector position as 'o' 
% Show actual end effector position as 'x'

function check_config
    %load trajectory
   % trajectory=trajectory(2:end,:); % Remove first line of values
    flag_prism = false;
    Fid = fopen('trajectory');
    no_links = fscanf(Fid, '%d', 1);
    for i=1:no_links
        trajectory(i,1) = fscanf(Fid, '%f',1);
        trajectory(i,2) = fscanf(Fid, '%f',1);
    end
    
    load angles
    %load arm
    Fid = fopen('arm');
    no_links = fscanf(Fid, '%d', 1);
    Lambda  = fscanf(Fid, '%f', 1);
    for i=1:no_links
        link_length(i) = fscanf(Fid, '%f',1);
        joint_var(i) = fscanf(Fid, '%f', 1);
        link_type(i) = fscanf(Fid, '%s',1);
                    
    end
            
    [m,n]=size(angles);
    
    dtheta=diff(angles);
    %link_length=arm(2:n+1,1);
    
    x_actual=zeros(m,2);
    error=zeros(m,1);
    mag_dtheta=zeros(m-1);
    for i=1:m-1
        mag_dtheta(i)=norm(dtheta(i,:));
    end
    
    
    for i=1:m
        x_actual(i,:)=[0,0];
        for j = 1:no_links
            if link_type(j) == 'P'
                flag_prism = true;
                joint_var = angles(i,j);
                link_length(j) = joint_var;
                angles(i,j) = 0;

            end
        end
        Sum=cumsum(angles(i,:))';
%       x_actual(i,1)=x_actual(i,1)+link_length'*cos(Sum);
%       x_actual(i,2)=x_actual(i,2)+link_length'*sin(Sum);
        x_actual(i,1)=x_actual(i,1)+link_length*cos(Sum);
        x_actual(i,2)=x_actual(i,2)+link_length*sin(Sum);
        %x_error=trajectory(i,:)-x_actual(i,:)';
        %error(i)=x_error*x_error';
    end
    
    figure(3)
        th=linspace(-pi/2,pi/2,1000);
        for i = 1:m
            text(trajectory(i,1)+0.015,trajectory(i,2),string(i))
             
        end
        hold on
        plot(trajectory(:,1),trajectory(:,2),'ob')
        hold on
        plot(x_actual(:,1),x_actual(:,2),'+r')
        hold on 
        if flag_prism == false
            R=ones(1,n)*link_length';
            plot(R*cos(th),R*sin(th),'g')
            legend('x_{desired}','x_{actual}','Workspace Boundary')
            hold on 
        else
            legend('x_{desired}','x_{actual}')
        end
        axis([min(trajectory(:,1))-0.5,max(trajectory(:,1))+0.5...
            min(trajectory(:,2))-0.5,max(trajectory(:,2))+0.5])
        title('Desired and Actual Position Trajectories')
        xlabel('x_1')
        ylabel('x_2')
        
end