classdef FlipManager < handle
  properties
    flip_x_axis
    flip_y_axis
    delta, omega;
    N;
    method;
  end


  methods

    function this= FlipManager(method, delta, N)

      this.delta = delta;
      this.omega = sqrt(9.81/0.3); % sqrt(g/h)
      this.N = N;
      x_axis = 1;
      y_axis = 2;
      this.method = method;
      this.flip_x_axis = Flipm( this.delta, this.omega, method, this.N, x_axis);
      this.flip_y_axis = Flipm( this.delta, this.omega, method, this.N, y_axis);
    end

    function cycle( this, numOfSteps )
      if nargin < 2
        numOfSteps = 200;
      end
      stateDim = 5;
      data = zeros( numOfSteps , stateDim*2+this.N*2);

      for i = 1:numOfSteps
        this.flip_x_axis.step();
        this.flip_y_axis.step();
        data(i, 1:stateDim)= this.flip_x_axis.q';
        data(i, stateDim+1:stateDim*2)= this.flip_y_axis.q';
        data(i, stateDim*2+1:stateDim*2+this.N)= this.flip_x_axis.z_dot';
        data(i, stateDim*2+this.N+1:stateDim*2+this.N*2)= this.flip_y_axis.z_dot';
      end

      fprintf(" computation ended \n")
      save('data.mat', 'data');
    end

    function plotSystemEvolution( this )

      dataStruct = load('data.mat');
      data = dataStruct.data;
      stateDim = (size(data,2) - this.N*2)/2;

      figure('Name','FLIPM System Evolution','pos',[10 10 1800 1150]);


      axisLimits_x = [-0.01 0.49 -0.13 0.13 ];           
      axisLimits_y = [-0.45 0.45 -0.3 0.3 ];           

      subplot(2,2,1); grid on;
      title('Side view');

        
      subplot(2,2,2); grid on;         
      title('Front view');

      for  i= 1:size(data,1)
        q_x = data(i,1:stateDim);
        z_dot_x= data(stateDim*2+1);
        subplot(2,2,1);

        drawings_x = Flipm.staticDraw( q_x, z_dot_x, axisLimits_x); 

        q_y = data(i,stateDim+1:stateDim*2);
        z_dot_y= data(stateDim*2+this.N+1);

        subplot(2,2,2);
        drawings_y = Flipm.staticDraw( q_y, z_dot_y,axisLimits_y);
        subplot(2,2,[3 4]);

        % state  q =  [ z, c1, c1_dot, c2, c2_dot ]
        zx = data(i,1);
        zy = data(i,stateDim+1);
        c1x = data(i,2);
        c1y = data(i,stateDim+2);
        c2x = data(i,4);
        c2y = data(i,stateDim+4);
        z_dot_x = data(i,stateDim*2+1:stateDim*2+this.N);
        z_dot_y = data(i,stateDim*2+this.N+1:stateDim*2+this.N*2);

        plot(data(1:i,1),data(1:i,stateDim+1), 'LineWidth', 1.0);
        hold on;
        plot(data(1:i,2),data(1:i,stateDim+2), 'LineWidth', 0.7);
        plot(data(1:i,4),data(1:i,stateDim+4), 'LineWidth', 0.4);
        plot(zx + this.delta*tril(ones(this.N,this.N))*z_dot_x', zy + ...
          this.delta*tril(ones(this.N,this.N))*z_dot_y', 'g', 'LineWidth', 1.0);
        axis equal;
        axis([-0.1,0.8,-0.2,0.2])
        legend('ZMP', 'm_1', 'm_2', 'predicted ZMP')
        grid on; xlabel('x [m]'), ylabel('y [m]');
        if ( this.method == "approximate")
            title('Approximate inverse');
        else
            title('Exact inverse');
        end

        hold off;
        pause( 0.0001);

        if i< size(data,1)
          Flipm.deleteDrawings( drawings_x);
          Flipm.deleteDrawings( drawings_y);
        end
      end
    end
  end
end
