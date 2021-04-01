%% Simulation of flow in 2D CPS disc
% https://www.mathworks.com/help/pde/ug/solve-poissons-equation-on-a-unit-disk.html

% Create the PDE with the geometry 
model = createpde();
geometryFromEdges(model,@circleg);

% Plot the geometry 
figure 
pdegplot(model,'EdgeLabels','on'); 
axis equal

% Zero Dirichlet B.C. on edges
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);

% Specify Coefficients 
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',1);

% Create a mesh with max element size 0.1
hmax = 0.1;
generateMesh(model,'Hmax',hmax);
figure
pdemesh(model); 
axis equal

% Solve and plot PDE 
results = solvepde(model);
u = results.NodalSolution;
pdeplot(model,'XYData',u)
title('Numerical Solution');
xlabel('x')
ylabel('y')

% Compare with analytical solution and plot error
p = model.Mesh.Nodes;
exact = (1 - p(1,:).^2 - p(2,:).^2)/4;
pdeplot(model,'XYData',u - exact')
title('Error');
xlabel('x')
ylabel('y')

% Refine mesh 
hmax = 0.1;
error = [];
err = 1;
while err > 5e-7 % run until error <= 5e-7
    generateMesh(model,'Hmax',hmax); % refine mesh
    results = solvepde(model);
    u = results.NodalSolution;
    p = model.Mesh.Nodes;
    exact = (1 - p(1,:).^2 - p(2,:).^2)/4; 
    err = norm(u - exact',inf); % compare with exact solution
    error = [error err]; % keep history of err
    hmax = hmax/2;
end

% Plot final mesh and solution
figure
pdemesh(model); 
axis equal
figure
pdeplot(model,'XYData',u)
title('Numerical Solution');
xlabel('x')
ylabel('y')




