% Create a phase space plot of the action potential V and restitution
% variable W.
% V0, W0: Initial values
% t_end and delta_t: Integration time and time step
% Call by, e.g., phaseplot( .3, 0., 50, 0.05 )
function phaseplot( V0, W0, t_end, delta_t )

clf();
% Create phase portrait
subplot(1, 2, 1);
portrait();

% Create array of discrete time values
n_step = ceil(t_end / delta_t);
t = delta_t * ( 0 : n_step );
t(end) = t_end;

% Create arrays for associated V and W values, set initial values
V = 0*t; W = 0*t;
V(1) = V0; W(1) = W0;

% Integrate up to t_end using Euler steps
for n = 2:numel(t)


end

% Plot solution in phase space
subplot(1, 2, 1);
plot( V, W, 'g*' )
% ... and vs. time
subplot(1, 2, 2);
plot(t, V, 'r');
hold on
plot(t, W, 'k');
legend('V', 'W');
end

function portrait()
% Phase space positions to plot
[ V_mesh, W_mesh ] = meshgrid(0 : .1 : 1, 0 : 0.2 : 2.5);
% Arrows with Vdot, Wdot
[ V_dot, W_dot ] = alpa(V_mesh, W_mesh);
quiver(V_mesh, W_mesh, V_dot, W_dot, 1.1);
hold on
% set plot limits
xlim( [-.1, 1.1] );
ylim( [-.1, 2.6] );
% set labels
xlabel('V')
ylabel('W')

end
