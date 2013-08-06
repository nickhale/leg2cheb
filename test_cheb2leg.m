
% Make a toy problem:
N = 1000;
rng(1);
c_cheb = randn(N+1, 1);

% Compute using asymptotic-based algorithm:
c_leg = cheb2leg(c_cheb);
% Compute using recurrence relation:
c_leg2 = cheb2leg(c_cheb, 0);
% Compute using Piessen's algorithm:
% c_leg2 = cheb2leg_piessens(c_cheb);

% The error:
err = norm(c_leg - c_leg2, inf)

