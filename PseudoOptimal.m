classdef PseudoOptimal < handle
    %PSEUDOOPTIMAL pseudo optimal control for simulation, parameter estimation, and optimal control
    %   dependency: chebfun, casadi, ChebTest
    
    properties
        npts;   % number of nodes
        nS;     % dim of states
        nU;     % dim of control
        nO;     % dim of output
        sGuess; % initial guess of states
        Pv0;    % initial value of variable parameters
        getP;   % function handle to get full parameter
        P;      % actually used parameter
        S;      % casadi fcn
    end
    
    properties (Hidden)
       orth;   % orthogonal nodes, weights, differential matrix
       guess;  % guess of initial states, states and parameter, or states and control
               % also time for nodes
       scale;  % time scale
       x;      % casadivariable for variable at nodes
       
       
       %Dynfcn; % Dynfcn:   dynamics function handle for vectorized discretized states and control
               %           format: Dynfcn( Xc,Uc,D,nS,scale,P )
               %                   Xc:    state at nodes
               %                   Uc:    control at nodes
               %                   D:     differential matrix
               %                   nS:    dim of states
               %                   scale: time scale
               %                   P:     parameter
       %Outfcn; % function handle of difference between output and target
               %           format: Outfcn( Xc )
               %           output nOutxN matrix(row for out, col for time)
        
       %target; % target of output, NxnO(row for time, col for var)
       %Pmin;   % lower bound of variable parameter
       %Pmax;   % upper bound of variable parameter
       %outerr; % acceptable difference between output and target
    end
    
    methods
        % function for simulation
        function [Xc,Xsim] = simulate(obj,Dynfcn,time,u,init,maxiter,displev)
            % simulate use collocation
            % requires:
            %   Dynfcn:   dynamics function handle for vectorized discretized states and control
            %             format: Dynfcn( Xc,Uc,D,nS,scale,P )
            %                     Xc:    state at nodes
            %                     Uc:    control at nodes
            %                     D:     differential matrix
            %                     nS:    dim of states
            %                     scale: time scale
            %                     P:     parameter
            %   time:     simulation time
            %   u:        input, corresponding to time(col for control, row for time)
            %   init:     initial condition for states
            %   maxiter:  maximum iteration for optimization
            %   displev:  disp level for optimization
            % output:
            %   Xc:       discretized states(col for state, row for time)
            %   Xsim:     simulated states corresponding to time(col for state, row for time)
            type = 'Simulation';
            import casadi.*;
            disp([type,': Initialize optimization problem ......']);
            init=init(:);% reshape to col vec
            time=time(:).';% reshape to row vec
            fullx=@(x)[init;x];% from xc to get full x with init cond
            [obj.orth.xx,obj.orth.ww,obj.orth.vv]=chebpts(obj.npts);% cheb pnts and weights
            obj.orth.D=ChebTest.getDifferentiationMatrix(obj.orth.xx,obj.orth.vv);% differential matrix
            obj.guess.tSpan=[time(1),time(end)];% get time span
            obj.guess.time=chebpts(obj.npts,obj.guess.tSpan).';% get time nodes
            ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end).');% initial of states at nodes
            obj.guess.state=reshape(ls.',numel(ls),1);% get initial state vec
            obj.scale=(obj.guess.tSpan(2)-obj.guess.tSpan(1))/2;% time scale
            Uc=interp1(time(:),u,obj.guess.time(:));% control at nodes
            Uc=reshape(Uc.',numel(Uc),1);
            eqn = @(x)Dynfcn(fullx(x),Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.P);% dynamics constraint
            
            obj.x=SX.sym('x',(obj.npts-1)*obj.nS,1);% casadi variable for state at node
            nlp=struct('x',obj.x,'f',dot(eqn(obj.x),eqn(obj.x)));% objective = minimize constraints
            opts.ipopt.max_iter=maxiter;% set maximum iteration number
            opts.ipopt.print_level=displev;% set display level
            obj.S=nlpsol('S','ipopt',nlp,opts);% optimization problem
            disp([type,': Initialize over, solving ......']);
            r=obj.S('x0',obj.guess.state);% solve optimization problem
            xopt=full(r.x);% get numerical result
            
            Xc = fullx(xopt);
            Xc = reshape(Xc,obj.nS,[]);
            Xc = Xc.';
            
            Xsim=ChebTest.barycentricInterpolate(time(:),Xc,obj.guess.time(:),obj.orth.vv(:));
        end
        
        % function for simulation
        function [Xc,Xsim] = simulatePre(obj,Dynfcn,time,u,init,maxiter,displev)
            % simulate use collocation
            % requires:
            %   Dynfcn:   dynamics function handle for vectorized discretized states and control
            %             format: Dynfcn( Xc,Uc,D,nS,scale,P )
            %                     Xc:    state at nodes
            %                     Uc:    control at nodes
            %                     D:     differential matrix
            %                     nS:    dim of states
            %                     scale: time scale
            %                     P:     parameter
            %   time:     simulation time
            %   u:        input, corresponding to time(col for control, row for time)
            %   init:     initial condition for states
            %   maxiter:  maximum iteration for optimization
            %   displev:  disp level for optimization
            % output:
            %   Xc:       discretized states(col for state, row for time)
            %   Xsim:     simulated states corresponding to time(col for state, row for time)
            type = 'Simulation';
            import casadi.*;
            disp([type,': Initialize optimization problem ......']);
            init=init(:);% reshape to col vec
            time=time(:).';% reshape to row vec
            fullx=@(x)[init;x];% from xc to get full x with init cond
            [obj.orth.xx,obj.orth.ww,obj.orth.vv]=chebpts(obj.npts);% cheb pnts and weights
            obj.orth.D=ChebTest.getDifferentiationMatrix(obj.orth.xx,obj.orth.vv);% differential matrix
            obj.guess.tSpan=[time(1),time(end)];% get time span
            obj.guess.time=chebpts(obj.npts,obj.guess.tSpan).';% get time nodes
            ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end).');% initial of states at nodes
            obj.guess.state=reshape(ls.',numel(ls),1);% get initial state vec
            obj.scale=(obj.guess.tSpan(2)-obj.guess.tSpan(1))/2;% time scale
            Uc=interp1(time(:),u,obj.guess.time(:));% control at nodes
            Uc=reshape(Uc.',numel(Uc),1);
            eqn = @(x)Dynfcn(fullx(x),Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.P);% dynamics constraint
            
            %obj.x=MX.sym('x',(obj.npts-1)*obj.nS,1);% casadi variable for state at node
            %nlp=struct('x',obj.x,'f',dot(eqn(obj.x),eqn(obj.x)));% objective = minimize constraints
            %opts.ipopt.max_iter=maxiter;% set maximum iteration number
            %opts.ipopt.print_level=displev;% set display level
            %S=nlpsol('S','ipopt',nlp,opts);% optimization problem
            disp([type,': Initialize over, solving ......']);
            r=obj.S('x0',obj.guess.state);% solve optimization problem
            xopt=full(r.x);% get numerical result
            
            Xc = fullx(xopt);
            Xc = reshape(Xc,obj.nS,[]);
            Xc = Xc.';
            
            Xsim=ChebTest.barycentricInterpolate(time(:),Xc,obj.guess.time(:),obj.orth.vv(:));
        end
        
        % function for parameter estimation
        function [Pvid,Xsim] = estimateParam(obj,Dynfcn,Outfcn,target,time,u,init,Pmin,Pmax,outerr,maxiter,displev)
            % parameter estimation use collocation
            % input:
            %   Dynfcn:   dynamics function handle for vectorized discretized states and control
            %             format: Dynfcn( Xc,Uc,D,nS,scale,P )
            %                     Xc:    state at nodes
            %                     Uc:    control at nodes
            %                     D:     differential matrix
            %                     nS:    dim of states
            %                     scale: time scale
            %                     P:     parameter
            %   Outfcn:   function handle output function
            %             format: Outfcn( Xc )
            %             output nOutxN matrix(row for out, col for time)
            %   target:   target of output, NxnO(row for time, col for var)
            %   time:     simulation time
            %   u:        input, corresponding to time(col for control, row for time)
            %   init:     initial condition for states
            %   outerr:   acceptable output error range
            %   Pmin:     lower bound of variable parameter
            %   Pmax:     upper bound of variable parameter
            %   maxiter:  maximum iteration for optimization
            %   displev:  disp level for optimization
            % output:
            %   Pvid:     discretized states(col for state, row for time)
            %   Xsim:     simulated states corresponding to time(col for state, row for time)
            type = 'Parameter Estimation';
            import casadi.*;
            disp([type,': Initialize optimization problem ......']);
            init=init(:);% reshape to col vec
            time=time(:).';% reshape to row vec
            fullx=@(x)[init;x];% from xc to get full x with init cond
            [obj.orth.xx,obj.orth.ww,obj.orth.vv]=chebpts(obj.npts);% cheb pnts and weights
            obj.orth.D=ChebTest.getDifferentiationMatrix(obj.orth.xx,obj.orth.vv);% differential matrix
            obj.guess.tSpan=[time(1),time(end)];% get time span
            obj.guess.time=chebpts(obj.npts,obj.guess.tSpan).';% get time nodes
            Targc=interp1(time,target,obj.guess.time(:));
            Targc=Targc.';% target at node
            ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end).');% initial of states at nodes
            obj.guess.state=reshape(ls.',numel(ls),1);% get initial state vec
            obj.guess.state=[obj.guess.state;obj.Pv0(:)];% get initial state and parameter vec
            obj.scale=(obj.guess.tSpan(2)-obj.guess.tSpan(1))/2;% time scale
            Uc=interp1(time(:),u,obj.guess.time(:));% control at nodes
            Uc=reshape(Uc.',numel(Uc),1);
            ind=obj.nS*(obj.npts-1);% index seperation of state and parameter
            %controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
            %    Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.getP(x(ind+1:end)));% dynamics constraint
            controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
                Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.getP(x(ind+1:end)));% dynamics constraint
            
            cost = @(x)dot(obj.orth.ww(:),1000*1/2*diag(x.'*x))*obj.scale;% cost for id
            %cost = @(x)dot(obj.orth.ww,sum(abs(x),1))*obj.scale;% 1-norm cost for id
            %cost = @(x)max(sum(abs(x),1));% inf-norm cost for id
            %cost = @(x)dot(obj.orth.ww,1/2*sum(x.^2,1))*obj.scale;% 2-norm cost for id
            nP = length(obj.Pv0);% dim of variable parameter
            
            controlieq = @(x)eye(nP)*x(ind+1:end);% parameter bound
            %controlieq2 = @(x)reshape((Outfcn(reshape(fullx(x(1:ind)),obj.nS,obj.npts))-Targc)*diag(1./obj.orth.ww),...
            %   obj.nO*obj.npts,1)/obj.scale;% output defect
            controlieq2 = @(x)reshape(Outfcn(reshape(fullx(x(1:ind)),obj.nS,obj.npts))-Targc,...
                obj.nO*obj.npts,1);% output defect
            
            obj.x=SX.sym('x',(obj.npts-1)*obj.nS+nP,1);% casadi variable for state at node and parameter
            nlp=struct('x',obj.x,'f',cost(Outfcn(reshape(fullx(obj.x(1:ind)),obj.nS,obj.npts))-Targc),...% objective = minimize difference between output and target
               'g',[controlieq(obj.x);controleqn(obj.x);controlieq2(obj.x)]);% constraints = parameter bound, dynamics, output defect
            
            opts.ipopt.max_iter=maxiter;% set maximum iteration number
            opts.ipopt.print_level=displev;% set display level
            obj.S=nlpsol('S','ipopt',nlp,opts);% optimization problem
            disp([type,': Initialize over, solving ......']);
            r=obj.S('x0',obj.guess.state,'lbg',[Pmin(:);-1e-5*ones(obj.npts*obj.nS,1);-outerr*ones(obj.nO*obj.npts,1)],...
               'ubg',[Pmax(:);1e-5*ones(obj.npts*obj.nS,1);outerr*ones(obj.nO*obj.npts,1)]);% solve optimization problem
            xopt=full(r.x);% get numerical result
            
            Xc = fullx(xopt(1:ind));
            Xc = reshape(Xc,obj.nS,[]);
            Xc = Xc.';
            
            Xsim=ChebTest.barycentricInterpolate(time(:),Xc,obj.guess.time(:),obj.orth.vv(:));
            Pvid = xopt(ind+1:end);
        end
        
        % function for parameter estimation
        function [Pvid,Xsim] = estimateParamInfn(obj,Dynfcn,Outfcn,target,time,u,init,Pmin,Pmax,outerr,maxiter,displev)
            % parameter estimation use collocation
            % input:
            %   Dynfcn:   dynamics function handle for vectorized discretized states and control
            %             format: Dynfcn( Xc,Uc,D,nS,scale,P )
            %                     Xc:    state at nodes
            %                     Uc:    control at nodes
            %                     D:     differential matrix
            %                     nS:    dim of states
            %                     scale: time scale
            %                     P:     parameter
            %   Outfcn:   function handle output function
            %             format: Outfcn( Xc )
            %             output nOutxN matrix(row for out, col for time)
            %   target:   target of output, NxnO(row for time, col for var)
            %   time:     simulation time
            %   u:        input, corresponding to time(col for control, row for time)
            %   init:     initial condition for states
            %   outerr:   acceptable output error range
            %   Pmin:     lower bound of variable parameter
            %   Pmax:     upper bound of variable parameter
            %   maxiter:  maximum iteration for optimization
            %   displev:  disp level for optimization
            % output:
            %   Pvid:     discretized states(col for state, row for time)
            %   Xsim:     simulated states corresponding to time(col for state, row for time)
            type = 'Parameter Estimation';
            import casadi.*;
            disp([type,': Initialize optimization problem ......']);
            init=init(:);% reshape to col vec
            time=time(:).';% reshape to row vec
            fullx=@(x)[init;x];% from xc to get full x with init cond
            [obj.orth.xx,obj.orth.ww,obj.orth.vv]=chebpts(obj.npts);% cheb pnts and weights
            obj.orth.D=ChebTest.getDifferentiationMatrix(obj.orth.xx,obj.orth.vv);% differential matrix
            obj.guess.tSpan=[time(1),time(end)];% get time span
            obj.guess.time=chebpts(obj.npts,obj.guess.tSpan).';% get time nodes
            Targc=interp1(time,target,obj.guess.time(:));
            Targc=Targc.';% target at node
            ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end).');% initial of states at nodes
            obj.guess.state=reshape(ls.',numel(ls),1);% get initial state vec
            obj.guess.state=[obj.guess.state;obj.Pv0(:)];% get initial state and parameter vec
            obj.scale=(obj.guess.tSpan(2)-obj.guess.tSpan(1))/2;% time scale
            Uc=interp1(time(:),u,obj.guess.time(:));% control at nodes
            Uc=reshape(Uc.',numel(Uc),1);
            ind=obj.nS*(obj.npts-1);% index seperation of state and parameter
            %controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
            %    Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.getP(x(ind+1:end)));% dynamics constraint
            controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
                Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.getP(x(ind+1:end)));% dynamics constraint
            
            %cost = @(x)dot(obj.orth.ww(:),1000*1/2*diag(x.'*x))*obj.scale;% cost for id
            %cost = @(x)dot(obj.orth.ww,sum(abs(x),1))*obj.scale;% 1-norm cost for id
            cost = @(x)norm(sum(abs(x),1),'inf');% inf-norm cost for id
            %cost = @(x)dot(obj.orth.ww,1/2*sum(x.^2,1))*obj.scale;% 2-norm cost for id
            nP = length(obj.Pv0);% dim of variable parameter
            
            controlieq = @(x)eye(nP)*x(ind+1:end);% parameter bound
            %controlieq2 = @(x)reshape((Outfcn(reshape(fullx(x(1:ind)),obj.nS,obj.npts))-Targc)*diag(1./obj.orth.ww),...
            %   obj.nO*obj.npts,1)/obj.scale;% output defect
            controlieq2 = @(x)reshape(Outfcn(reshape(fullx(x(1:ind)),obj.nS,obj.npts))-Targc,...
                obj.nO*obj.npts,1);% output defect
            
            obj.x=SX.sym('x',(obj.npts-1)*obj.nS+nP,1);% casadi variable for state at node and parameter
            nlp=struct('x',obj.x,'f',cost(Outfcn(reshape(fullx(obj.x(1:ind)),obj.nS,obj.npts))-Targc),...% objective = minimize difference between output and target
               'g',[controlieq(obj.x);controleqn(obj.x);controlieq2(obj.x)]);% constraints = parameter bound, dynamics, output defect
            
            opts.ipopt.max_iter=maxiter;% set maximum iteration number
            opts.ipopt.print_level=displev;% set display level
            obj.S=nlpsol('S','ipopt',nlp,opts);% optimization problem
            disp([type,': Initialize over, solving ......']);
            r=obj.S('x0',obj.guess.state,'lbg',[Pmin(:);-1e-5*ones(obj.npts*obj.nS,1);-outerr*ones(obj.nO*obj.npts,1)],...
               'ubg',[Pmax(:);1e-5*ones(obj.npts*obj.nS,1);outerr*ones(obj.nO*obj.npts,1)]);% solve optimization problem
            xopt=full(r.x);% get numerical result
            
            Xc = fullx(xopt(1:ind));
            Xc = reshape(Xc,obj.nS,[]);
            Xc = Xc.';
            
            Xsim=ChebTest.barycentricInterpolate(time(:),Xc,obj.guess.time(:),obj.orth.vv(:));
            Pvid = xopt(ind+1:end);
        end
        
        % function for parameter estimation
        function [Pvid,Xsim] = estimateParamPre(obj,Dynfcn,Outfcn,target,time,u,init,Pmin,Pmax,outerr,Pv0)
            % parameter estimation use collocation
            % input:
            %   Dynfcn:   dynamics function handle for vectorized discretized states and control
            %             format: Dynfcn( Xc,Uc,D,nS,scale,P )
            %                     Xc:    state at nodes
            %                     Uc:    control at nodes
            %                     D:     differential matrix
            %                     nS:    dim of states
            %                     scale: time scale
            %                     P:     parameter
            %   Outfcn:   function handle output function
            %             format: Outfcn( Xc )
            %             output nOutxN matrix(row for out, col for time)
            %   target:   target of output, NxnO(row for time, col for var)
            %   time:     simulation time
            %   u:        input, corresponding to time(col for control, row for time)
            %   init:     initial condition for states
            %   outerr:   acceptable output error range
            %   Pmin:     lower bound of variable parameter
            %   Pmax:     upper bound of variable parameter
            %   maxiter:  maximum iteration for optimization
            %   displev:  disp level for optimization
            % output:
            %   Pvid:     discretized states(col for state, row for time)
            %   Xsim:     simulated states corresponding to time(col for state, row for time)
            type = 'Parameter Estimation';
            import casadi.*;
            disp([type,': Initialize optimization problem ......']);
            init=init(:);% reshape to col vec
            time=time(:).';% reshape to row vec
            fullx=@(x)[init;x];% from xc to get full x with init cond
            [obj.orth.xx,obj.orth.ww,obj.orth.vv]=chebpts(obj.npts);% cheb pnts and weights
            obj.orth.D=ChebTest.getDifferentiationMatrix(obj.orth.xx,obj.orth.vv);% differential matrix
            obj.guess.tSpan=[time(1),time(end)];% get time span
            obj.guess.time=chebpts(obj.npts,obj.guess.tSpan).';% get time nodes
            Targc=interp1(time,target,obj.guess.time(:));
            Targc=Targc.';% target at node
            ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end).');% initial of states at nodes
            obj.guess.state=reshape(ls.',numel(ls),1);% get initial state vec
            %obj.guess.state=[obj.guess.state;obj.Pv0(:)];% get initial state and parameter vec
            obj.guess.state=[obj.guess.state;Pv0(:)];% get initial state and parameter vec
            obj.scale=(obj.guess.tSpan(2)-obj.guess.tSpan(1))/2;% time scale
            Uc=interp1(time(:),u,obj.guess.time(:));% control at nodes
            Uc=reshape(Uc.',numel(Uc),1);
            ind=obj.nS*(obj.npts-1);% index seperation of state and parameter
            %controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
            %    Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.getP(x(ind+1:end)));% dynamics constraint
            controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
                Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.getP(x(ind+1:end)));% dynamics constraint
            
            %cost = @(x)dot(obj.orth.ww(:),1000*1/2*diag(x.'*x))*obj.scale;% cost for id
            cost = @(x)dot(obj.orth.ww,sum(abs(x),1))*obj.scale;% 1-norm cost for id
            %cost = @(x)dot(obj.orth.ww,1/2*sum(x.^2,1))*obj.scale;% 2-norm cost for id
            nP = length(obj.Pv0);% dim of variable parameter
            
            controlieq = @(x)eye(nP)*x(ind+1:end);% parameter bound
            %controlieq2 = @(x)reshape((Outfcn(reshape(fullx(x(1:ind)),obj.nS,obj.npts))-Targc)*diag(1./obj.orth.ww),...
            %   obj.nO*obj.npts,1)/obj.scale;% output defect
            controlieq2 = @(x)reshape(Outfcn(reshape(fullx(x(1:ind)),obj.nS,obj.npts))-Targc,...
                obj.nO*obj.npts,1);% output defect
            
            %obj.x=MX.sym('x',(obj.npts-1)*obj.nS+nP,1);% casadi variable for state at node and parameter
            %nlp=struct('x',obj.x,'f',cost(Outfcn(reshape(fullx(obj.x(1:ind)),obj.nS,obj.npts))-Targc),...% objective = minimize difference between output and target
            %   'g',[controlieq(obj.x);controleqn(obj.x);controlieq2(obj.x)]);% constraints = parameter bound, dynamics, output defect
            
            %opts.ipopt.max_iter=maxiter;% set maximum iteration number
            %opts.ipopt.print_level=displev;% set display level
            %S=nlpsol('S','ipopt',nlp,opts);% optimization problem
            disp([type,': Initialize over, solving ......']);
            r=obj.S('x0',obj.guess.state,'lbg',[Pmin(:);-1e-5*ones(obj.npts*obj.nS,1);-outerr*ones(obj.nO*obj.npts,1)],...
               'ubg',[Pmax(:);1e-5*ones(obj.npts*obj.nS,1);outerr*ones(obj.nO*obj.npts,1)]);% solve optimization problem
            xopt=full(r.x);% get numerical result
            
            Xc = fullx(xopt(1:ind));
            Xc = reshape(Xc,obj.nS,[]);
            Xc = Xc.';
            
            Xsim=ChebTest.barycentricInterpolate(time(:),Xc,obj.guess.time(:),obj.orth.vv(:));
            Pvid = xopt(ind+1:end);
        end
        
        % function for parameter estimation
        function [Pvid,Xsim] = estimateParamPreInit(obj,Dynfcn,Outfcn,target,time,u,init,Pmin,Pmax,outerr,sGuess,Pv0)
            % parameter estimation use collocation
            % input:
            %   Dynfcn:   dynamics function handle for vectorized discretized states and control
            %             format: Dynfcn( Xc,Uc,D,nS,scale,P )
            %                     Xc:    state at nodes
            %                     Uc:    control at nodes
            %                     D:     differential matrix
            %                     nS:    dim of states
            %                     scale: time scale
            %                     P:     parameter
            %   Outfcn:   function handle output function
            %             format: Outfcn( Xc )
            %             output nOutxN matrix(row for out, col for time)
            %   target:   target of output, NxnO(row for time, col for var)
            %   time:     simulation time
            %   u:        input, corresponding to time(col for control, row for time)
            %   init:     initial condition for states
            %   outerr:   acceptable output error range
            %   Pmin:     lower bound of variable parameter
            %   Pmax:     upper bound of variable parameter
            %   maxiter:  maximum iteration for optimization
            %   displev:  disp level for optimization
            % output:
            %   Pvid:     discretized states(col for state, row for time)
            %   Xsim:     simulated states corresponding to time(col for state, row for time)
            type = 'Parameter Estimation';
            import casadi.*;
            disp([type,': Initialize optimization problem ......']);
            init=init(:);% reshape to col vec
            time=time(:).';% reshape to row vec
            fullx=@(x)[init;x];% from xc to get full x with init cond
            [obj.orth.xx,obj.orth.ww,obj.orth.vv]=chebpts(obj.npts);% cheb pnts and weights
            obj.orth.D=ChebTest.getDifferentiationMatrix(obj.orth.xx,obj.orth.vv);% differential matrix
            obj.guess.tSpan=[time(1),time(end)];% get time span
            obj.guess.time=chebpts(obj.npts,obj.guess.tSpan).';% get time nodes
            Targc=interp1(time,target,obj.guess.time(:));
            Targc=Targc.';% target at node
            %ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end).');% initial of states at nodes
            ls = interp1(time(:),sGuess,obj.guess.time(2:end).');% initial of states at nodes
            obj.guess.state=reshape(ls.',numel(ls),1);% get initial state vec
            %obj.guess.state=[obj.guess.state;obj.Pv0(:)];% get initial state and parameter vec
            obj.guess.state=[obj.guess.state;Pv0(:)];% get initial state and parameter vec
            obj.scale=(obj.guess.tSpan(2)-obj.guess.tSpan(1))/2;% time scale
            Uc=interp1(time(:),u,obj.guess.time(:));% control at nodes
            Uc=reshape(Uc.',numel(Uc),1);
            ind=obj.nS*(obj.npts-1);% index seperation of state and parameter
            %controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
            %    Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.getP(x(ind+1:end)));% dynamics constraint
            controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
                Uc,obj.orth.D,obj.nS,obj.nU,obj.scale,obj.getP(x(ind+1:end)));% dynamics constraint
            
            %cost = @(x)dot(obj.orth.ww(:),1000*1/2*diag(x.'*x))*obj.scale;% cost for id
            cost = @(x)dot(obj.orth.ww,sum(abs(x),1))*obj.scale;% 1-norm cost for id
            %cost = @(x)dot(obj.orth.ww,1/2*sum(x.^2,1))*obj.scale;% 2-norm cost for id
            nP = length(obj.Pv0);% dim of variable parameter
            
            controlieq = @(x)eye(nP)*x(ind+1:end);% parameter bound
            %controlieq2 = @(x)reshape((Outfcn(reshape(fullx(x(1:ind)),obj.nS,obj.npts))-Targc)*diag(1./obj.orth.ww),...
            %   obj.nO*obj.npts,1)/obj.scale;% output defect
            controlieq2 = @(x)reshape(Outfcn(reshape(fullx(x(1:ind)),obj.nS,obj.npts))-Targc,...
                obj.nO*obj.npts,1);% output defect
            
            %obj.x=MX.sym('x',(obj.npts-1)*obj.nS+nP,1);% casadi variable for state at node and parameter
            %nlp=struct('x',obj.x,'f',cost(Outfcn(reshape(fullx(obj.x(1:ind)),obj.nS,obj.npts))-Targc),...% objective = minimize difference between output and target
            %   'g',[controlieq(obj.x);controleqn(obj.x);controlieq2(obj.x)]);% constraints = parameter bound, dynamics, output defect
            
            %opts.ipopt.max_iter=maxiter;% set maximum iteration number
            %opts.ipopt.print_level=displev;% set display level
            %S=nlpsol('S','ipopt',nlp,opts);% optimization problem
            disp([type,': Initialize over, solving ......']);
            %r=obj.S('x0',obj.guess.state,'lbg',[Pmin(:);-1e-5*ones(obj.npts*obj.nS,1);-outerr*ones(obj.nO*obj.npts,1)],...
            %   'ubg',[Pmax(:);1e-5*ones(obj.npts*obj.nS,1);outerr*ones(obj.nO*obj.npts,1)]);% solve optimization problem
            r=obj.S('x0',obj.guess.state,'lbg',[Pmin(:);-0*1e-5*ones(obj.npts*obj.nS,1);kron(ones(obj.npts,1),-outerr)],...
               'ubg',[Pmax(:);0*1e-5*ones(obj.npts*obj.nS,1);kron(ones(obj.npts,1),outerr)]);% solve optimization problem
            xopt=full(r.x);% get numerical result
            
            Xc = fullx(xopt(1:ind));
            Xc = reshape(Xc,obj.nS,[]);
            Xc = Xc.';
            
            Xsim=ChebTest.barycentricInterpolate(time(:),Xc,obj.guess.time(:),obj.orth.vv(:));
            Pvid = xopt(ind+1:end);
        end
        
        % function for optimal control
        function [Xc,Uc,Xopt,Uopt] = optimalControl(obj,Dynfcn,Costfcn,Confcn,conlb,conub,time,u,init,maxiter,displev)
            % optimal control use collocation
            % input:
            %   Dynfcn:   dynamics function handle for vectorized discretized states and control
            %             format: Dynfcn( Xc,Uc,D,nS,scale,P )
            %                     Xc:    state at nodes
            %                     Uc:    control at nodes
            %                     D:     differential matrix
            %                     nS:    dim of states
            %                     scale: time scale
            %                     P:     parameter
            %   Costfcn:  function handle cost function
            %             format: Costfcn( Xc,Uc,tc,D,nS,nU,w,scale )
            %   Confcn:   function handle of constraint function !!!: need modification and consideration...
            %             format: Confcn( Xc,Uc,D,scale )
            %   conlb:    lower bound of constraint function
            %   conub:    upper bound of constraint function
            %   target:   target of output, NxnO(row for time, col for var)
            %   time:     simulation time
            %   u:        input, corresponding to time(col for control, row for time)
            %   init:     initial condition for states
            %   maxiter:  maximum iteration for optimization
            %   displev:  disp level for optimization
            % output:
            %   Pvid:     discretized states(col for state, row for time)
            %   Xsim:     simulated states corresponding to time(col for state, row for time)
            type = 'Optimal Control';
            import casadi.*;
            disp([type,': Initialize optimization problem ......']);
            init=init(:);% reshape to col vec
            time=time(:).';% reshape to row vec
            fullx=@(x)[init;x];% from xc to get full x with init cond
            [obj.orth.xx,obj.orth.ww,obj.orth.vv]=chebpts(obj.npts);% cheb pnts and weights
            obj.orth.D=ChebTest.getDifferentiationMatrix(obj.orth.xx,obj.orth.vv);% differential matrix
            obj.guess.tSpan=[time(1),time(end)];% get time span
            obj.guess.time=chebpts(obj.npts,obj.guess.tSpan).';% get time nodes
            ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end).');% initial of states at nodes
            obj.guess.state=reshape(ls.',numel(ls),1);% get initial state vec
            Uc=interp1(time(:),u,obj.guess.time(:));% control at nodes
            Uc=reshape(Uc.',numel(Uc),1);
            obj.guess.state=[obj.guess.state;Uc];% get initial state and control vec
            obj.scale=(obj.guess.tSpan(2)-obj.guess.tSpan(1))/2;% time scale            
            ind=obj.nS*(obj.npts-1);% index seperation of state and control
            controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
                x(ind+1:end),obj.orth.D,obj.nS,obj.nU,obj.scale,obj.P);% dynamics constraint
            
            cost = @(x)Costfcn( fullx(x(1:ind)),x(ind+1:end),obj.guess.time,...
                                obj.orth.D,obj.nS,obj.nU,...
                                obj.orth.ww,obj.scale );% cost for optimal control
            controlieq = @(x)Confcn(fullx(x(1:ind)),x(ind+1:end),...
                             obj.orth.D,obj.scale);% constraintfcn bound
            
            obj.x=SX.sym('x',(obj.npts-1)*obj.nS+obj.npts*obj.nU,1);% casadi variable for state at node and parameter
            nlp=struct('x',obj.x,'f',cost(obj.x),...% objective = minimize cost fcn
                'g',[controlieq(obj.x);controleqn(obj.x)]);% constraints = extconstr, dynamics
            
            opts.ipopt.max_iter=maxiter;% set maximum iteration number
            opts.ipopt.print_level=displev;% set display level
            obj.S=nlpsol('S','ipopt',nlp,opts);% optimization problem
            disp([type,': Initialize over, solving ......']);
            r=obj.S('x0',obj.guess.state,'lbg',[kron(ones(obj.npts,1),conlb);zeros(obj.npts*obj.nS,1)],...
                'ubg',[kron(ones(obj.npts,1),conub);zeros(obj.npts*obj.nS,1)]);% solve optimization problem
            xopt=full(r.x);% get numerical result
            
            Xc = fullx(xopt(1:ind));
            Xc = reshape(Xc,obj.nS,[]);
            Xc = Xc.';
            Uc = xopt(ind+1:end);
            Uc = reshape(Uc,obj.nU,[]);
            Uc = Uc.';
            
            Xopt=ChebTest.barycentricInterpolate(time(:),Xc,obj.guess.time(:),obj.orth.vv(:));
            Uopt=ChebTest.barycentricInterpolate(time(:),Uc,obj.guess.time(:),obj.orth.vv(:));
        end
        
        % function for optimal control
        function [Xc,Uc,Xopt,Uopt] = optimalControlTerminal(obj,Dynfcn,Costfcn,Confcn,conlb,conub,time,u,init,terminal,maxiter,displev)
            % optimal control use collocation
            % input:
            %   Dynfcn:   dynamics function handle for vectorized discretized states and control
            %             format: Dynfcn( Xc,Uc,D,nS,scale,P )
            %                     Xc:    state at nodes
            %                     Uc:    control at nodes
            %                     D:     differential matrix
            %                     nS:    dim of states
            %                     scale: time scale
            %                     P:     parameter
            %   Costfcn:  function handle cost function
            %             format: Costfcn( Xc,Uc,tc,D,nS,nU,w,scale )
            %   Confcn:   function handle of constraint function !!!: need modification and consideration...
            %             format: Confcn( Xc,Uc,D,scale )
            %   conlb:    lower bound of constraint function
            %   conub:    upper bound of constraint function
            %   target:   target of output, NxnO(row for time, col for var)
            %   time:     simulation time
            %   u:        input, corresponding to time(col for control, row for time)
            %   init:     initial condition for states
            %   maxiter:  maximum iteration for optimization
            %   displev:  disp level for optimization
            % output:
            %   Pvid:     discretized states(col for state, row for time)
            %   Xsim:     simulated states corresponding to time(col for state, row for time)
            type = 'Optimal Control';
            import casadi.*;
            disp([type,': Initialize optimization problem ......']);
            init=init(:);% reshape to col vec
            terminal=terminal(:);
            time=time(:).';% reshape to row vec
            %fullx=@(x)[init;x];% from xc to get full x with init cond
            fullx=@(x)[init;x;terminal];% from xc to get full x with init cond
            [obj.orth.xx,obj.orth.ww,obj.orth.vv]=chebpts(obj.npts);% cheb pnts and weights
            obj.orth.D=ChebTest.getDifferentiationMatrix(obj.orth.xx,obj.orth.vv);% differential matrix
            obj.guess.tSpan=[time(1),time(end)];% get time span
            obj.guess.time=chebpts(obj.npts,obj.guess.tSpan).';% get time nodes
            %ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end).');% initial of states at nodes
            ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end-1).');% initial of states at nodes
            obj.guess.state=reshape(ls.',numel(ls),1);% get initial state vec
            Uc=interp1(time(:),u,obj.guess.time(:));% control at nodes
            Uc=reshape(Uc.',numel(Uc),1);
            obj.guess.state=[obj.guess.state;Uc];% get initial state and control vec
            obj.scale=(obj.guess.tSpan(2)-obj.guess.tSpan(1))/2;% time scale            
            %ind=obj.nS*(obj.npts-1);% index seperation of state and control
            ind=obj.nS*(obj.npts-2);% index seperation of state and control
            controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
                x(ind+1:end),obj.orth.D,obj.nS,obj.nU,obj.scale,obj.P);% dynamics constraint
            
            cost = @(x)Costfcn( fullx(x(1:ind)),x(ind+1:end),obj.guess.time,...
                                obj.orth.D,obj.nS,obj.nU,...
                                obj.orth.ww,obj.scale );% cost for optimal control
            controlieq = @(x)Confcn(fullx(x(1:ind)),x(ind+1:end),...
                             obj.orth.D,obj.scale);% constraintfcn bound
            
            %obj.x=SX.sym('x',(obj.npts-1)*obj.nS+obj.npts*obj.nU,1);% casadi variable for state at node and parameter
            obj.x=SX.sym('x',(obj.npts-2)*obj.nS+obj.npts*obj.nU,1);% casadi variable for state at node and parameter
            nlp=struct('x',obj.x,'f',cost(obj.x),...% objective = minimize cost fcn
                'g',[controlieq(obj.x);controleqn(obj.x)]);% constraints = extconstr, dynamics
            
            opts.ipopt.max_iter=maxiter;% set maximum iteration number
            opts.ipopt.print_level=displev;% set display level
            obj.S=nlpsol('S','ipopt',nlp,opts);% optimization problem
            disp([type,': Initialize over, solving ......']);
            r=obj.S('x0',obj.guess.state,'lbg',[kron(ones(obj.npts,1),conlb);zeros(obj.npts*obj.nS,1)],...
                'ubg',[kron(ones(obj.npts,1),conub);zeros(obj.npts*obj.nS,1)]);% solve optimization problem
            xopt=full(r.x);% get numerical result
            
            Xc = fullx(xopt(1:ind));
            Xc = reshape(Xc,obj.nS,[]);
            Xc = Xc.';
            Uc = xopt(ind+1:end);
            Uc = reshape(Uc,obj.nU,[]);
            Uc = Uc.';
            
            Xopt=ChebTest.barycentricInterpolate(time(:),Xc,obj.guess.time(:),obj.orth.vv(:));
            Uopt=ChebTest.barycentricInterpolate(time(:),Uc,obj.guess.time(:),obj.orth.vv(:));
        end
        
        % function for optimal control
        function [Xc,Uc,Xopt,Uopt] = optimalControlPre(obj,Dynfcn,Costfcn,Confcn,conlb,conub,time,u,init,maxiter,displev,sGuess)
            % optimal control use collocation
            % input:
            %   Dynfcn:   dynamics function handle for vectorized discretized states and control
            %             format: Dynfcn( Xc,Uc,D,nS,scale,P )
            %                     Xc:    state at nodes
            %                     Uc:    control at nodes
            %                     D:     differential matrix
            %                     nS:    dim of states
            %                     scale: time scale
            %                     P:     parameter
            %   Costfcn:  function handle cost function
            %             format: Costfcn( Xc,Uc,tc,D,nS,nU,w,scale )
            %   Confcn:   function handle of constraint function !!!: need modification and consideration...
            %             format: Confcn( Xc,Uc,D,scale )
            %   conlb:    lower bound of constraint function
            %   conub:    upper bound of constraint function
            %   target:   target of output, NxnO(row for time, col for var)
            %   time:     simulation time
            %   u:        input, corresponding to time(col for control, row for time)
            %   init:     initial condition for states
            %   maxiter:  maximum iteration for optimization
            %   displev:  disp level for optimization
            % output:
            %   Pvid:     discretized states(col for state, row for time)
            %   Xsim:     simulated states corresponding to time(col for state, row for time)
            type = 'Optimal Control';
            import casadi.*;
            disp([type,': Initialize optimization problem ......']);
            init=init(:);% reshape to col vec
            time=time(:).';% reshape to row vec
            fullx=@(x)[init;x];% from xc to get full x with init cond
            [obj.orth.xx,obj.orth.ww,obj.orth.vv]=chebpts(obj.npts);% cheb pnts and weights
            obj.orth.D=ChebTest.getDifferentiationMatrix(obj.orth.xx,obj.orth.vv);% differential matrix
            obj.guess.tSpan=[time(1),time(end)];% get time span
            obj.guess.time=chebpts(obj.npts,obj.guess.tSpan).';% get time nodes
            %ls = interp1(time(:),obj.sGuess,obj.guess.time(2:end).');% initial of states at nodes
            ls = interp1(time(:),sGuess,obj.guess.time(2:end).');% initial of states at nodes
            obj.guess.state=reshape(ls.',numel(ls),1);% get initial state vec
            Uc=interp1(time(:),u,obj.guess.time(:));% control at nodes
            Uc=reshape(Uc.',numel(Uc),1);
            obj.guess.state=[obj.guess.state;Uc];% get initial state and control vec
            obj.scale=(obj.guess.tSpan(2)-obj.guess.tSpan(1))/2;% time scale            
            ind=obj.nS*(obj.npts-1);% index seperation of state and control
            controleqn = @(x)Dynfcn(fullx(x(1:ind)),...
                x(ind+1:end),obj.orth.D,obj.nS,obj.nU,obj.scale,obj.P);% dynamics constraint
            
            cost = @(x)Costfcn( fullx(x(1:ind)),x(ind+1:end),obj.guess.time,...
                                obj.orth.D,obj.nS,obj.nU,...
                                obj.orth.ww,obj.scale );% cost for optimal control
            controlieq = @(x)Confcn(fullx(x(1:ind)),x(ind+1:end),...
                             obj.orth.D,obj.scale);% constraintfcn bound
            
            obj.x=SX.sym('x',(obj.npts-1)*obj.nS+obj.npts*obj.nU,1);% casadi variable for state at node and parameter
            %nlp=struct('x',obj.x,'f',cost(obj.x),...% objective = minimize cost fcn
            %    'g',[controlieq(obj.x);controleqn(obj.x)]);% constraints = extconstr, dynamics
            
            opts.ipopt.max_iter=maxiter;% set maximum iteration number
            opts.ipopt.print_level=displev;% set display level
            %obj.S=nlpsol('S','ipopt',nlp,opts);% optimization problem
            disp([type,': Initialize over, solving ......']);
            r=obj.S('x0',obj.guess.state,'lbg',[kron(ones(obj.npts,1),conlb);zeros(obj.npts*obj.nS,1)],...
                'ubg',[kron(ones(obj.npts,1),conub);zeros(obj.npts*obj.nS,1)]);% solve optimization problem
            xopt=full(r.x);% get numerical result
            
            Xc = fullx(xopt(1:ind));
            Xc = reshape(Xc,obj.nS,[]);
            Xc = Xc.';
            Uc = xopt(ind+1:end);
            Uc = reshape(Uc,obj.nU,[]);
            Uc = Uc.';
            
            Xopt=ChebTest.barycentricInterpolate(time(:),Xc,obj.guess.time(:),obj.orth.vv(:));
            Uopt=ChebTest.barycentricInterpolate(time(:),Uc,obj.guess.time(:),obj.orth.vv(:));
        end
        
    end
    
end

