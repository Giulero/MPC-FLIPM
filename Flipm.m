classdef Flipm < handle
  properties
    %number of future steps of npc
    N;
    % state  q =  [ z, c1, c1_dot, c2, c2_dot ]
    q;
    z_dot; % desired zmp
    za; % measured zmp 
    zmpMax; zmpMin;
    constraints;
    axis; % used to choose the zmp constraint
    % matrices we need
    A, B, C, D, % system matrices
    % flip properties
    m1, m2, b, k, tau;
    d_max; % max distance between c1 and c2
    %mpc 
    %intrinsic stability coefficients
    delta, omega;
    % constraint matrices
    A_inq, b_inq;
    A_eq, b_eq;
    C_MPC;
  end

 methods

    function this = Flipm( delta, omega, method, N, axis)

      this.N = N;
      this.axis = axis;
      this.q= zeros(5,1);
      this.za = 0;

      this.m1 = 4; this.m2 = 0.5;
      this.b = 20; this.k = 1000;
      this.tau = 0.00000001;
      this.d_max = 0.1; % the max distance betwee c1 and c2
      this.delta = delta;
      this.omega = omega;

      % if we decrease the b the d_max must be increased

      this.z_dot = zeros(this.N, 1);

      this.constraints = [Constraint([0,0], [0.04,0.02], 20), ...
        Constraint([0.0,0.1], [0.04,0.02], 20), ...
        Constraint([0.1,-0.1], [0.04,0.02], 20), ...
        Constraint([0.2,0.1], [0.04,0.02], 20), ...
        Constraint([0.3,-0.1], [0.04,0.02], 20), ...
        Constraint([0.4,0.1], [0.04,0.02], 20), ...
        Constraint([0.5,0], [0.04,0.02], 1000)];

      this.computeZmpBounds();
      this.computeSystemMatrices( method);
      this.computeMpcMatrices();
    end

    function computeSystemMatrices(this, method)

      M = (this.m1+this.m2)/(this.m1*this.m2);

      if ( method == "approximate")


        alpha = 1/(this.m1*this.m2)*(this.tau*(this.k - M*this.b^2) + this.b);
        beta = this.k/(this.m1*this.m2)*(1 - M*this.b*this.tau);

        D_a = beta*this.tau/(this.m1*this.m2);

        a = this.omega^2/(D_a*this.m2);
        b = this.k/this.m2 + beta/D_a;
        c = this.b/this.m2 + alpha/D_a;

       A = [0 0 0 0 0;
             0 0 1 0 0;
             0 -this.k/this.m1 -this.b/this.m1 this.k/this.m1 this.b/this.m1;
             0 0 0 0 1;
            -a (a+b) c -b -c];
            
        B = [1 0 0 0 0]';
         
      else % exact inverse

        a = (this.b*M*this.omega^2)*(this.m1/this.b); % x des
        b = this.omega^2*(this.m1/this.b); % x_dot des
        c = -a; % zmp des
        d = this.k/this.m2;
        e = this.b/this.m2 + this.k/this.b;
        f = -d;
        g = -e;
        ctr = -b;
        A = [0 0 0 0 0;
             0 0 1 0 0;
             0 -this.k/this.m1 -this.b/this.m1 this.k/this.m1 this.b/this.m1;
             0 0 0 0 1;
             c (d+a) (e+b) f g];

        B = [1 0 0 0 ctr]';
           
      end
      % not the ZMP but the distance between the two masses
      C = [ 0 -1 0 1 0];
      D = 0;
           
      sys = ss(A,B,C,D);
      d_sys = c2d(sys, this.delta);
      this.A = d_sys.A;
      this.B = d_sys.B;
      this.C = d_sys.C;
      this.D = d_sys.D;
    end

    function computeMpcMatrices( this )
      state_dim = size(this.A,1);
      this.A_inq = zeros(4*this.N,this.N);
      this.b_inq = zeros(4*this.N,1);
      this.C_MPC = zeros(this.N, state_dim);
      % boundedness condition
      this.A_eq = ...
        (1-exp(-this.omega*this.delta))/ ...
        (this.omega*(1-exp(-this.omega*this.delta*this.N)))* ...
        exp(-this.omega*this.delta*(0:this.N-1));
    
      % fill the matrix for the "distance between the two masses" constraint
      for i=1:this.N
        for j=1:i
          this.A_inq(i, j) = this.C*this.A^(i-j)*this.B;
        end
        this.C_MPC(i,:) = this.C*this.A^(i);
      end
      
      % stack "distance" constraint
      this.A_inq(this.N+1:2*this.N, :) = -this.A_inq(1:this.N, :);
      % stack ZMP constraint
      this.A_inq(2*this.N+1:end, :) = [ 
      (this.delta*tril(ones(this.N, this.N)));
      -(this.delta*tril(ones(this.N, this.N)));
      ];
    end

    function step( this )
      this.computeZmpBounds();
      this.computeMpcInput();
      this.updateState();
      this.shiftConstraints();
    end

    function computeMpcInput( this )

      % Mpc iteration:   cost function, just for one axis
      % think about to modify in order to embed the dynamic of the system
      H = eye(this.N);

      this.b_inq(1:2*this.N) = [ 
      (this.d_max - this.C_MPC*this.q);
      (this.d_max + this.C_MPC*this.q)];
      z = this.q(1);
      % vector containing the zmp reference
      z_des = z + this.delta*tril(ones(this.N, this.N))* this.z_dot;


      this.b_inq(2*this.N+1:end) = [ 
          (this.zmpMax(1:this.N) - z);
         -(this.zmpMin(1:this.N) - z)
         ];
      % Compute the equality constraint
      this.b_eq = this.q(2) + this.q(3)/this.omega - z;

      options = optimset('Algorithm','interior-point-convex','Display','off');

      [z_dot,~,exitflag,~]= ...
        quadprog(H,[],this.A_inq,this.b_inq,this.A_eq,this.b_eq,[],[],[],options);

      if exitflag == -2
        fprintf("---SOLUTION NOT FOUND---");
      end

      this.z_dot = z_dot;
 
    end

    function updateState(this )
     this.q = this.A*this.q + this.B*this.z_dot(1);

     this.za = this.q(2)  + 1/(this.m1*this.omega^2)* ...
     (this.b*(this.q(3) - this.q(5)) + this.k*(this.q(2) -this.q(4)));
    end

    function computeZmpBounds(this)
      %% computes the ZmpBoundaries. this.axis specifies if x
      % boundaries or y boundaries.
      this.zmpMax = [];
      this.zmpMin = [];

      for i = 1:size(this.constraints,2)
        this.zmpMax = [this.zmpMax; (this.constraints(i).center(this.axis) + ...
          this.constraints(i).size(this.axis))*ones(this.constraints(i).duration,1)];
        this.zmpMin = [this.zmpMin; (this.constraints(i).center(this.axis) - ...
          this.constraints(i).size(this.axis))*ones(this.constraints(i).duration,1)];
      end
    end

    function shiftConstraints(this)
      % Shift first constraint back in time, if not active anymore, delete it
      this.constraints(1).duration = this.constraints(1).duration - 1;
      if this.constraints(1).duration == 0
        this.constraints(1) = [];
      end
    end

  end
  methods(Static = true)
    function  drawings = staticDraw( q, z_dot, axisLimits  )
      h = 0.4;
      d = Drawer();
      orange =[0.88,0.45,0.02];
      green=[0.66,0.88,0.02];
      turquoise=[0.02,0.88,0.88];
      blue=[0.02,0.45,0.88];
      red = [ 0.88, 0.02, 0.02];
      black = [ 0, 0, 0];


      c1= [ q(2); 0 ];
      c2= [ q(4); 0 ];

      hold on;

      drawings = [
      d.drawCircle2D( c1(1), c1(2), 0.002, red);
      d.drawCircle2D( c2(1), c2(2), 0.001, black);
      ];

      axis(axisLimits);
      legend('c1', 'c2');
    end
    function deleteDrawings( drawings)
      for i = 1:size(drawings,1)
        delete(drawings(i,1));
      end
    end
  end
end

