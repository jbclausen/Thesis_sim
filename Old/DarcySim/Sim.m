%% Jens Clausen - Darcy Law Simulation 
% let x-y plane be the plane on which the disk is. The z axis represents
% the time taken for the flow to reach that point on the disc. 

% Define constants 
mu = 8.9e-4;
gamma = 72.86e-3;
alpha = 71;
r_p = 10e-6;
r_c = 2.54e-3;

% Set up grid 
x=linspace(-3,3,40);
y=linspace(-3,3,40);
[x,y]=meshgrid(x,y);
z = zeros(40,40);

% Calculate time to the edge of the disc
sum = 0;
for i = 1:(size(x,1)/2)
    for j = 1:(size(y,1)/2)
        r = sqrt((x(i)*x(i))+(y(j)*y(j)));
        % z(i,j) = r; % test - store radius in z
        
        
    end 
end

% Draw mesh 
mesh(x,y,z)