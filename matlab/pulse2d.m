function pulse2d

% Use same grid spacing delta_x along x and y, but different
% #intervals
M_x = 200;
L_x = 300;
M_y = 200;

t_end = 5000.;
plot_interval = 100.;

% Diffusivity:
eta = .3;

delta_x = L_x / M_x
L_y = delta_x * M_y

xl = (0:M_x) * delta_x;
yl = (0:M_y) * delta_x;

[x, y] = meshgrid(xl, yl);
 
delta_t = 0.2 * delta_x^2 / eta
n_step = ceil(t_end / delta_t);
t = delta_t * ( 0 : n_step );
t(end) = t_end;

V = 0*x;
W = V;

% spiral:


% pulse:


show(xl, yl, V, W, 0.)
next_plot = plot_interval;

for n = 2:numel(t)

    % Euler step
    
    % Neumann at y_min/y_max (spiral):

    % Dirichlet at y_min/y_max (pulse):
       %V(1, :) = 0; V(end, :) = 0.;

    % Neumann at x_min/x_max:

    
    if t(n) >= next_plot
        next_plot = next_plot + plot_interval;
        show( xl, yl, V, W, t(n) );
    end
end
end

function show(x, y, V, W, t)
subplot(1, 2, 1)
image(x, y, 250*V);
tit = ['V @ t = ', num2str(t)];
title(tit);
xlabel('x');
ylabel('y');
subplot(1, 2, 2)
image(x, y, 50*W);
title('W')
pause(.2)
end
