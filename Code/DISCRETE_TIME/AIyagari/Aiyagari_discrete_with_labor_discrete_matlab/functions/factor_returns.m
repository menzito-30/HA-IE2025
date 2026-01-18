function [R, W] = factor_returns(KN, par, grid)

r = par.alpha     * KN^(par.alpha-1) - par.delta;
R = 1 + r;
W = (1-par.alpha) * KN^par.alpha ;