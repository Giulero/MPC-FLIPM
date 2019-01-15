%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC_based control for FLIP                    %
% Course: Underactuated Robotics                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% change flag to execute one method or the other
method = "approximate"; % approximate inverse
method = "exact";

delta = 0.02; % time step
N = 50; % preview window
fM = FlipManager('exact', delta, N);
%also on console
fM.cycle(); %main computation
fM.plotSystemEvolution(); %graphical front end
