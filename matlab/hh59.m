function hh59

delta_t = 0.05;
t_end = 10.;
nStep = floor( t_end / delta_t );

clf

V = zeros(1, nStep);
n = V; m = V; h =  V; t = V;

V(1) = 7;
t(1) = 0;

n(1) = alphan(0) / ( alphan(0) + betan(0)) ;
% ...

times = 5.1;
ntstep = times*nStep;

for i = 1:nStep
    
% ...
    
t(i+1) = t(i) + delta_t;

end

subplot(3, 1, 1);
plot(t,V, '-'); legend('V(t)', 'Location', 'NorthEast');
subplot(3, 1, 2);
plot(t, Ik(n, V), t, Ina(m, h, V), t, Il(V) );
legend('I Ka^+', 'I Na^+', 'I leak', 'Location', 'NorthEast');
subplot(3, 1, 3);
plot(t, n.^4, t, m.^3, t, h );
legend('Gate Ka^+ n', 'Gate Na^+ m^3', 'Gate Na^+ h', 'Location', 'NorthEast');

hold off
end

function ret = VG(V,n,m,h)
C = 1;
ret = -1/C * ( Ik(n,V) + Il(V) + Ina(m,h,V) );
end

function ret = Il(V)
gl = 0.3;
Vl = 10;
ret = gl*(V-Vl);
end

function ret = Ik(n,V)
gk = 36;
Vk = -12;
ret = gk * n.^4 .* (V-Vk);
end

function ret = Ina(m,h,V)
gna = 120;
Vna = 115;
ret = gna * m.^3 .* h .*(V-Vna);
end

function ret = mG(m,V)
ret = alpham(V)*(1-m) - betam(V)*m;
end

function ret = hG(h,V)
ret = alphah(V)*(1-h) - betah(V)*h;
end

function ret = nG(n,V)
ret = alphan(V)*(1-n) - betan(V)*n;
end

function ret = alphan(V)
ret = 0.01*(V - 10) ./ (1 - exp( (10-V)/10) );
end

function ret = alphah(V)
ret = 0.07*exp(-V/20);
end

function ret = alpham(V)
ret = 0.1*(V - 25) ./ (1 - exp( (25-V)/10) );
end

function ret = betan(V)
faca = 0.125;
ret = faca * exp(-V/80);
end

function ret = betam(V)
ret = 4 * exp(-V/18);
end

function ret = betah(V)
ret = 1 ./ ( 1 + exp( (30 - V)/10) );
end
