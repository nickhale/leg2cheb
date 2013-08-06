
% Make a toy problem:
N = 1000;
rng(1);
c_leg = randn(N, 1);

% Compute using asymptotic-based algorithm:
c_cheb = leg2cheb(c_leg);
% Compute using recurrence relation:
c_cheb2 = leg2cheb(c_leg, 0);
% Compute using Piessen's algorithm:
% c_cheb2 = leg2cheb_piessens(c_leg);

% The error:
err = norm(c_cheb - c_cheb2, inf)
