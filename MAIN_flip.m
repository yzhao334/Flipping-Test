% MAIN_flip.m
% remember to add casadi into path

clc; clear;
% setting parameters
drone_params;
% Trajectory Parameters:
% global N z0 zF state_target control_target;
N=10; % N should be even
duration = 1/200*N*2*10;
% duration = 1/200*N;

% Initial State:
z0 = [0 ;0; -0.46; 0; 0; 0; 0; 0; 0; 0; 0; 0];  % initial state

% zF = [0 ;0; -0.46; 0; 0; 0; 0; 0; 0; 0; 0; 0];  % final state
zF = [0 ;0; -0.46; 2*pi; 0; 0; 0; 0; 0; 0; 0; 0];  % final state

state_target=nan(12,N*2);
temptime=linspace(0,duration,N*2);
for i=1:6
    state_target(i,1:N/2)=z0(i)*ones(1,N/2);
    state_target(i,N/2+1:3/2*N)=linspace(z0(i),zF(i),N);
    state_target(i,3/2*N+1:2*N)=zF(i)*ones(1,N/2);
end
for i=7:9
    state_target(i,:)=gradient(state_target(i-6,:))./gradient(temptime);
end
% state_target(10,:)=gradient(state_target(6,:))./gradient(temptime);
% state_target(11,:)=gradient(state_target(5,:))./gradient(temptime);
% state_target(12,:)=gradient(state_target(4,:))./gradient(temptime);
state_target(10,:)=gradient(state_target(5,:))./gradient(temptime);
state_target(11,:)=gradient(state_target(4,:))./gradient(temptime);
state_target(12,:)=gradient(state_target(6,:))./gradient(temptime);


ct=quad.Ct*quad.rho*quad.A*quad.r^2;
control_target=ones(4,N*2)*sqrt(quad.M*quad.g/ct/4);
% control_target=inv_dyn(state_target,temptime,quad);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( droneDynamics(x,u,quad) );
problem.func.pathObj = @(t,x,u)( objective(x,u) );

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Set up bounds on state and control                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = duration;
problem.bounds.finalTime.upp = duration;

problem.bounds.initialState.low = z0;
problem.bounds.initialState.upp = z0;
problem.bounds.finalState.low = zF;
problem.bounds.finalState.upp = zF;


problem.bounds.control.low = sqrt(10*quadEDT.motorcommandToW2_gain)*[1;1;1;1]*0;
problem.bounds.control.upp = sqrt(quadEDT.motors_max*quadEDT.motorcommandToW2_gain)*[1;1;1;1];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Initialize trajectory with guess                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Car travels at a speed of one, and drives in a straight line from start
% to finish point.

% problem.guess.time = [0,duration];
% problem.guess.state = [problem.bounds.initialState.low, problem.bounds.finalState.low];
% problem.guess.control = [zeros(4,1),zeros(4,1)];

% problem.guess.time = linspace(0,duration,N*2);
% problem.guess.state=state_target;
% problem.guess.control=control_target;

problem.guess.time = [0,duration/4,duration/4*3,duration];
problem.guess.state=[z0,z0,zF,zF];
problem.guess.control=control_target(:,1)*ones(1,4);

% problem.guess.time = [0,duration];
% problem.guess.state=[z0,z0];
% problem.guess.control=control_target(:,1)*ones(1,2);
    

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e5);


problem.options.method = 'trapezoid';
% problem.options.method = 'hermiteSimpson';
% problem.options.method = 'rungeKutta';
% problem.options.method = 'chebyshev';

% problem.options.method = 'chebyshev';
% problem.options.chebyshev.nColPts = N*2;
problem.options.trapezoid.nGrid = N*2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
sim.u0=[-0.66708 0 0 0].';

neo.tspan=[0,0.95];
neo.dt=1e-3;
neo.time=0:neo.dt:neo.tspan(end);
neo.init=state_target(:,1);
neo.terminal=state_target(:,end);
neo.t_tar=linspace(neo.tspan(1),neo.tspan(end),size(state_target,2));
neo.target=state_target.';
neo.dyn=@(Xc,Uc,D,nS,nU,scale,P)dyntest(Xc,Uc,problem.func.dynamics,D,scale);
neo.cost=@(Xc,Uc,tc,D,nS,nU,ww,scale)costTest(Xc,Uc,...
    neo.target,neo.t_tar,tc,ww,scale,sim.u0);
neo.constr=@(Xc,Uc,D,scale)Uc;
neo.conlb=problem.bounds.control.low;
% neo.conub=problem.bounds.control.upp;
% neo.conlb=problem.bounds.control.upp*0.03;
neo.conub=problem.bounds.control.upp*0.85;
% neo.conlb=-sqrt(quadEDT.motorcommandToW2_gain*500)*ones(4,1)*0.85;
% neo.conub=+sqrt(quadEDT.motorcommandToW2_gain*500)*ones(4,1)*0.85;

ct=quad.Ct*quad.rho*quad.A*quad.r^2;
cQ=quad.Cq*quad.rho*quad.A*quad.r^3;
d=sqrt(2)/2*quad.d;
M =...
[  -ct,   -ct,   -ct,   -ct;...
  ct*d, -ct*d, -ct*d,  ct*d;...
  ct*d,  ct*d, -ct*d, -ct*d;...
   -cQ,    cQ,   -cQ,    cQ];
% M_inv_neg=-inv(M);

sim.u_init=[(sim.u0+[-0.63;0;2e-2;0])*ones(1,34),...            
            (sim.u0+[-0.63;0;-2e-2;0])*ones(1,34),...            
            (sim.u0+[-0.63;0;0;0])*ones(1,40),...
            (sim.u0+[0.2;0;0;0])*ones(1,42)];
sim.u_init=M\sim.u_init;
sim.u_init=sqrt(sim.u_init);
sim.t_temp=linspace(neo.tspan(1),neo.tspan(end),size(sim.u_init,2));

sim.u_init=sqrt(quad.M*quad.g/ct/4)*ones(size(sim.u_init));
neo.u_init=interp1(sim.t_temp(:),sim.u_init.',neo.time(:));

%%
ps=PseudoOptimal;
ps.npts=20;
ps.nS=12;
ps.nU=4;
ps.sGuess=interp1(neo.t_tar(:),neo.target,neo.time(:));
%%
[Xc,Uc,Xopt,Uopt] = ps.optimalControlTerminal(neo.dyn,neo.cost,neo.constr,...
    neo.conlb,neo.conub,neo.time,neo.u_init,...
    neo.init,neo.terminal,1000,5);

%%
% ps.npts=30;

%%
ps.sGuess=Xopt;
neo.u_init=Uopt;

[Xc,Uc,Xopt,Uopt] = ps.optimalControlTerminal(neo.dyn,neo.cost,neo.constr,...
    neo.conlb,neo.conub,neo.time,neo.u_init,...
    neo.init,neo.terminal,1000,5);

% soln = trajOpt(problem);

% get results
% t=linspace(soln.grid.time(1),soln.grid.time(end),2*N);
% z=soln.interp.state(t);
% u=soln.interp.control(t);
% t=soln.grid.time;
% z=soln.grid.state;
% u=soln.grid.control;
%% simulate
% load staticGain;
% K=static.K;
sim.dt=1/200;
sim.soln.time=0;
sim.soln.state=neo.init;
sim.T=neo.time(end)/sim.dt;
sim.dyn=@(t,x,u,p)droneDynamics(x,u,p);
sim.soln.u=[];

[t,y]=ode23tb(@(t,x)sim.dyn(t,x,(interp1(neo.time,Uopt,t)).',quad),neo.tspan,neo.init);
sim.time=neo.tspan(1):sim.dt:neo.tspan(end);
sim.soln.state=(interp1(t(:),y,sim.time(:))).';

% for k=0:sim.T-1
%     sim.soln.time=[sim.soln.time,(k+1)*sim.dt];
%     % use opt control, Euler diff
%     %temp_s=(interp1(neo.time,Xopt,k*sim.dt)).';
%     %u_fb=sqrt(0*abs(quadEDT.motorcommandToW2_gain*(K*(temp_s-sim.soln.state(:,end)))));
%     %u_fb=0*([5*ones(4,6),10*ones(4,6)]*(temp_s-sim.soln.state(:,end)));
%     u_fb=0;
%     sim.soln.u=[sim.soln.u,(interp1(neo.time,Uopt,k*sim.dt)).'+u_fb];
%     temp.y=sim.soln.state(:,end)+sim.dt*feval(sim.dyn,0,sim.soln.state(:,end),sim.soln.u(:,end),quad);
%     sim.soln.state=[sim.soln.state,temp.y(:,end)];
% end
% clear k temp;
%% plot
sim.anim_freq=100; % refresh freq for animation
quad_anim(sim.soln.state,quad,sim.anim_freq);% animation
