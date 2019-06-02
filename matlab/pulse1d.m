% One-dimensional monodomain computation
function pulse1d

% Build x-Array that divides [0, L] into M intervals
M = 100;
L = 200;

t_end = 8000.;
plot_interval = 50.;

% Diffusivity:
eta = .3;

delta_x = L / M
x = (0:M) * delta_x;

% time step from CFL stability condition
delta_t = 0.1 * delta_x^2 / eta
n_step = ceil(t_end / delta_t);
t = delta_t * ( 0 : n_step );
t(end) = t_end;

% Initial values for V, W
V = 0*x;
W = V;

show(x, V, W, 0.)
next_plot = plot_interval;

% variables for CV and BCL detection
previous_V = 0.;
threshold = .1;
upTime = -10.;

for n = 2:numel(t)
    dV_diff = ...
    [ dV_reac, dW_reac ] = alpa(V, W);
    % Euler step. Apply time scale of 12.9 ms from Aliev-Panfilov model to
    % reaction terms.
    dt = t(n) - t(n-1);
    V = ...
    W = ...

    
    jetzt = t(n);
    if jetzt >= next_plot
        show(x, V, W, jetzt);
        next_plot = next_plot + plot_interval;
    end
    
    % Measure CV and BCL by monitoring one position
    current_V = V(10);
    if( (current_V > threshold) && (previous_V <= threshold) )
        % upstroke detected
        CV = L / ( jetzt - upTime )
        BCL = jetzt - upTime
        upTime = jetzt;
    end
    if( (current_V < threshold) && (previous_V >= threshold) )
        % decay detected
        APD = jetzt - upTime
     end
    previous_V = current_V;
end
end

function show(x, v, w, t)
hold off
plot(x, v, 'r-'); hold on;
plot(x, .5*w, 'b-'); hold on;
xlabel('x'); ylabel('f');
ylim([0, 1.1]);
title( ['t = ', num2str(t) ] );
pause(0.2)
end
