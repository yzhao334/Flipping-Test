# Flipping-Test
Test for solving drone flipping problem using numerical optimal control (specific, pseudo optimal control)

Fast example: run MAIN_flip.m (fast result)

Long example: run MAIN_flip_long.m (result include stabilizing stage after flipping)

dependent packages: chebfun-master(v5.7.0), casadi-matlabR2014b-v3.2.3(linux x64).

File chain: MAIN_flip -- costTest, MAIN_flip_long -- costTest_long

public files: 
              
              (plotting) draw_drone, draw_update,quad_anim
              (dynamics) droneDynamics dyntest drone_params
              (optimal control) PseudoOptimal ChebTest

Acknowlegment:

               Trajectory Optimization (package) - Matthew Kelly
               "Multirotor aerial vehicles." Mahony, Robert, Vijay Kumar, and Peter Corke. IEEE Robotics and Automation magazine 20.32 (2012)               
               Rolling Spider software package for Education (MIT toolbox)
               ME290J 2016 spring lectures, UC Berkeley, Francesco Borrelli

===

License:    BSD  (https://github.com/yzhao334/Flipping-Test/blob/master/LICENSE)
