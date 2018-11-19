delta = 0.03;
omega = sqrt(9.81/0.3); % sqrt(g/h)
N = 50;

contr_x = Controller(delta, omega, N, 1);
contr_y = Controller(delta, omega, N, 2);

cdesx = 0; cdesy = 0;
c1x = 0; c1xd = 0;
c1y = 0; c1yd = 0;
c2x = 0; c2xd = 0;
c2y = 0; c2yd = 0;
zx = 0; zy = 0;
z_meas_x = 0; z_meas_y = 0;

for i = 1:500
    
    fprintf("iteration n.%i\n", i)
    clf
    contr_x.computeMpc();
    contr_y.computeMpc();
    
    cdesx(end+1) = contr_x.c1_des;
    cdesy(end+1) = contr_y.c1_des;
    
    c1x(end+1) = contr_x.c1;
    c1y(end+1) = contr_y.c1;
    
    c2x(end+1) = contr_x.c2;
    c2y(end+1) = contr_y.c2;
    
    c1xd(end+1) = contr_x.c1_dot;
    c1yd(end+1) = contr_y.c1_dot;
    
    zx(end+1) = contr_x.z;
    zy(end+1) = contr_y.z;
    z_meas_x(end+1) = contr_x.za;
    z_meas_y(end+1) = contr_y.za;   
    
    hold on

    plot(zx, zy, 'r')
    plot(zx(end-1) + delta*tril(ones(N,N))*contr_x.z_dot, zy(end-1) + delta*tril(ones(N,N))*contr_y.z_dot, 'g')
    plot(c1x, c1y , 'm')
    plot(cdesx, cdesy , 'y');
    plot(contr_x.c1_des, contr_y.c1_des, '-s')
    plot(contr_x.c1, contr_y.c1, 'o')
    plot(c2x, c2y, 'b' );
    plot(z_meas_x, z_meas_y) 
    plot(c1x(end) + c1xd(end)/omega, c1y(end) + c1yd(end)/omega, '*') 
    axis equal
    axis([-0.1,0.8,-0.2,0.2])
    legend('zmp_{act_{des}}', 'zmp_{fut}', 'com1_{traj}', 'com_{traj_{des}}', 'com1_{des}', 'com1', 'com2', 'zmp_{meas}', 'CP')
%     printFootSteps(contr_x.constraints, i, N);
    drawnow
    if contr_x.constraints(end).duration == 0
        break
    end
end
