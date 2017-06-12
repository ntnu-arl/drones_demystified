%%  RRT grow example
%   This example is created as the matlab implementation of the classic RRT
%   implementation from Lavalle
%   (http://msl.cs.uiuc.edu/~lavalle/sub/rrt.py)

clear;
xdim = 640;
ydim = 480;
winsize = [xdim ydim];
epsilon_ = 7;
num_nodes = 5000;

nodes(1,1:2) = [xdim/2 ydim/2];

close all; figure; hold on; xlim([0 xdim]); ylim([0 ydim]); axis equal; grid on; box on;
for i = 1:num_nodes
    rand_point_ = [rand*xdim rand*ydim];
    nn = nodes(1,1:2);
    [n_nodes,~] = size(nodes);
    for j = 1:n_nodes
        if euclidean_dist(nodes(j,:),rand_point_) < euclidean_dist(nn,rand_point_)
            nn = nodes(j,:);
        end
    end
    newnode_ = step_from_to(nn,rand_point_,epsilon_);
    nodes = [newnode_; nodes];
    line([nn(1) newnode_(1)],[nn(2) newnode_(2)]);  drawnow;
 end
