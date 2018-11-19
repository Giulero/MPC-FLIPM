classdef Controller < dynamicprops
    % Alternative class for the controller
    % We try to combine the model of the LIP and the flip
    % in order to embed the input from mpc e constrain also the
    % distence between c1 and c2
    % Also split x and y controller
    
    properties
        delta, omega, N,
        c1, c1_dot, c1_ddot, % the actual CoM
        c1_des, c1_des_dot, % the desider CoM trajectory, from LIP.
        c2, c2_dot, % the actual small mass
        d3, d3_dot, % the distance betwee c1 and c2
        z; z_dot; % desired zmp
        za; % measured zmp 
        zmpMax; zmpMin;
        constraints;
        axis; % used to choose the zmp constraint
        % matrices we need
        A, B, C, D, % system matrices
        A_inq, b_inq, A_eq, b_eq, C_MPC % constraint matrices
        % flip properties
        m1, m2, b, k, tau;
        d_max; % max distance between c1 and c2
        H, F1, F2; % matrices for the MPC
        
        
    end
    
    methods
        function this = Controller(delta, omega, N, axis)
            %% initialization
            this.delta = delta;
            this.omega = omega;
            this.N = N;
            this.axis = axis;
            
            this.c1_des = 0; this.c1_des_dot = 0;
            this.c1 = 0.0; this.c1_dot = 0.0; this.c1_ddot = 0.0;
            this.c2 = 0.0; this.c2_dot = 0.0;
            this.z = 0; this.d3 = this.c2-this.c1; this.d3_dot = this.c2_dot-this.c1_dot;
            this.za = 0;
            
            % flip init --- TODO: FIND THE RIGHT VALUES!!!
            this.m1 = 4; this.m2 = 0.5;
            this.b = 20; this.k = 1000; this.tau = 0.00000001; 
            this.d_max = 0.1; % the max distance between c1 and c2
            
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
            this.computeMatrices();
            this.computeMPCmatrices();
        end
        
        function computeZmpBounds(this)
            %% computes the ZmpBoundaries. this.axis specifies if x
            % boundaries or y boundaries.
            this.zmpMax = [];
            this.zmpMin = [];
            
            for i = 1:size(this.constraints,2)
                this.zmpMax = [this.zmpMax; (this.constraints(i).center(this.axis) + this.constraints(i).size(this.axis))*ones(this.constraints(i).duration,1)];
                this.zmpMin = [this.zmpMin; (this.constraints(i).center(this.axis) - this.constraints(i).size(this.axis))*ones(this.constraints(i).duration,1)];
            end
        end
        
        function computeMpc(this)
            this.computeZmpBounds();
            %% Mpc iteration
            % cost function, just for one axis
            H = eye(this.N); % think about to modify in order to embed the dynamic of the system

            this.b_inq(1:2*this.N) = [(this.d_max - this.C_MPC*[this.c1_des, this.c1_des_dot, this.z, this.c1, this.c1_dot, this.c2, this.c2_dot, this.za]');
                                      (this.d_max + this.C_MPC*[this.c1_des, this.c1_des_dot, this.z, this.c1, this.c1_dot, this.c2, this.c2_dot, this.za]')];
                                 
            z = this.z;

            z_des = this.z + this.delta*tril(ones(this.N, this.N))*this.z_dot; % vector containing the zmp reference
                        
            x0 = [this.c1_des, this.c1_des_dot, this.z, this.c1, this.c1_dot, this.c2, this.c2_dot, this.za];
            
            f = x0*this.F1 + z_des'*this.F2;
            
            this.b_inq(2*this.N+1:end) = [(this.zmpMax(1:this.N) - z);
                                         -(this.zmpMin(1:this.N) - z)];
            
            % Compute the equality constraint
            this.b_eq = this.c1_des + this.c1_des_dot/this.omega - z;
            
            options = optimset('Algorithm','interior-point-convex','Display','off');
            
            [this.z_dot,~,exitflag,~]= quadprog(H,[],this.A_inq,this.b_inq,this.A_eq,this.b_eq,[],[],[],options);
            
            if exitflag == -2
                fprintf("---SOLUTION NOT FOUND---");
            end
            

            this.updateState()
            
            % Shift first constraint back in time, if not active anymore, delete it
            this.constraints(1).duration = this.constraints(1).duration - 1;
            if this.constraints(1).duration == 0
                this.constraints(1) = [];
            end
        end
        
        function computeMPCmatrices(this)
            % note that size(this.A)(1) is the number of the variables into
            % the state.
            
            this.A_inq = zeros(4*this.N,this.N); %
            this.b_inq = zeros(4*this.N,1);
            this.C_MPC = zeros(this.N,size(this.A,1));
            this.A_eq = (1-exp(-this.omega*this.delta))/(this.omega*...
                        (1-exp(-this.omega*this.delta*this.N)))*...
                        exp(-this.omega*this.delta*(0:this.N-1));

            for i=1:this.N
                for j=1:i
                    this.A_inq(i, j) = this.C*this.A^(i-j)*this.B;
                end
                this.C_MPC(i,:) = this.C*this.A^(i);
            end
            
            this.A_inq(this.N+1:2*this.N, :) = -this.A_inq(1:this.N, :);
            this.A_inq(2*this.N+1:end, :) = [(this.delta*tril(ones(this.N, this.N)));...
                                            -(this.delta*tril(ones(this.N, this.N)))];
            
            Q = 0*eye(size(this.A,1));
            R = 0.1;
            R_t = 0.01*eye(this.N);
            %% used for minimizinf the state 7
%             Q_t = eye(size(this.A,1)*this.N);
%             S_t = zeros(size(this.A,1)*this.N, this.N);
%             for i=1:this.N
%                 for j=1:i
%                     S_t((size(this.A,1)*(i-1)+1):(size(this.A,1)*(i-1)+size(this.A,1)),j) = this.A^(i-j)*this.B;
%                 end
%             end
%             T_t = zeros(size(this.A,1)*this.N, size(this.A,1));
%             for i=1:this.N
%                 T_t((size(this.A,1)*(i-1)+1):(size(this.A,1)*(i-1)+size(this.A,1)),:) = this.A^i;
%             end
            %% for minimizing the output state 8
            S_t = zeros(this.N, this.N);
            Q_t = 1*eye(this.N);
            C = [0 0 0 0 0 0 0 1];
            for i=1:this.N
                for j=1:i
                    S_t(i,j) = C*this.A^(i-j)*this.B;
                end
            end
            T_t = zeros(this.N, size(this.A,1));
            for i=1:this.N
                T_t(i,:) = C*this.A^i;
            end
            this.H = 2*(R_t+(S_t)'*Q_t*(S_t));
            this.H = (this.H+this.H')/2; % if Hessian matrix for the optimization probl. is not symmetric
            this.F1 = 2*((T_t)'*Q_t*(S_t));
            this.F2 = -2*(Q_t*S_t);
            Y = 2*(Q+T_t'*Q_t*T_t);
        end
               
        function computeMatrices(this)
            %% the state here is composed by (c, cd, z)des and (c1, c1d, c2, c2d, z_a)
            
            M = (this.m1+this.m2)/(this.m1*this.m2);
            
            A_2 = [        0         1;
                   -this.k*M -this.b*M];

            B_2 = [        0;
                   1/this.m2];

            alpha = 1/(this.m1*this.m2)*(this.tau*(this.k - M*this.b^2) + this.b);
            beta = this.k/(this.m1*this.m2)*(1 - M*this.b*this.tau);
%             fprintf("beta:%f \n", beta)

            C_a = [beta * this.m1 alpha * this.m2];
            D_a = beta*this.tau/(this.m1*this.m2);
%             fprintf("D_a:%f \n", D_a)
%             A = [ 0 1 0 0 0 0 0;
%                  this.omega^2 0 -this.omega^2 0 0 0 0;
%                  0 0 0 0 0 0 0;
%                  0 0 0 0 1 0 0;
%                  0 0 0 -this.k/this.m1 -this.b/this.m1 this.k/this.m1 this.b/this.m1;
%                  0 0 0 0 0 0 1;
%                  (this.omega^2/D_a) 0 -(this.omega^2/D_a) (this.k/this.m2 + beta/D_a*this.m2) (this.b + alpha/D_a*this.m2) -(this.k/this.m2 + beta/D_a*this.m2) -(this.b + alpha/D_a*this.m2)];
            A = [0 1 0 0 0 0 0;
                 this.omega^2 0 -this.omega^2 0 0 0 0;
                 0 0 0 0 0 0 0;
                 0 0 0 0 1 0 0;
                 0 0 0 -this.k/this.m1 -this.b/this.m1 this.k/this.m1 this.b/this.m1;
                 0 0 0 0 0 0 1;
                 (this.omega^2/(D_a*this.m2)) 0 -(this.omega^2/(D_a*this.m2)) (this.k/this.m2 + beta/D_a) (this.b/this.m2 + alpha/D_a) -(this.k/this.m2 + beta/D_a) -(this.b/this.m2 + alpha/D_a)];
            
            B = [0 0 1 0 0 0 0]';
            C = [0 0 0 -1 0 1 0];
            D = 0;
            
            sys = ss(A,B,C,D);
            d_sys = c2d(sys, this.delta);
            % don't get confused: the following are the matrices of the
            % discrete system!
            this.A = zeros(8,8);
            this.A(1:7,1:7) = d_sys.A;
            % increase the state with the zmp_meas at time k+1
            this.A(end,:) = this.A(4,:) - 1/(this.omega^2)*(this.k/this.m1*(this.A(6,:)-this.A(4,:))+...
                              this.b/this.m1*(this.A(7,:)-this.A(5,:)));
%             disp(this.k/this.m1*(this.A(6,6))+this.b/this.m1*(this.A(7,6)));
%             disp(this.A);
            this.B = zeros(8,1);
            this.B(1:7,:) = d_sys.B;
            this.C = zeros(1,8);
            this.C(:,1:7) = d_sys.C;
            this.D = d_sys.D;
            
            Co = ctrb(this.A, this.B);

            
            if rank(Co) == 8
                fprintf("Controllable system.\n");
            else
                fprintf("Non controllable\n");
                disp(Co);
                fprintf("Rank= %i\n", rank(Co));
%                 [~,a]=rref(Co)
            end
            
            for i=1:8
               C = zeros(1,8);
               C(i) = 1;
               fprintf("ctrb matrix-element %i has rank= %i \n", i,  checkOutputControllability(this.A, this.B, C));
            end

            
        end
        
        function updateState(this)
            %% the state here is composed by (c, cd, z)des and (c1, c1d, c2, c2d, z_a)
            
            x = this.A*[this.c1_des, this.c1_des_dot, this.z, this.c1, this.c1_dot, this.c2, this.c2_dot, this.za]' + this.B*this.z_dot(1);
            this.c1_des = x(1);
            this.c1_des_dot = x(2);
            this.z = x(3);
            this.c1 = x(4);
            this.c1_dot = x(5);
            this.c2 = x(6);
            this.c2_dot = x(7);
            this.za = x(8);
            
        end
        
    end
end

