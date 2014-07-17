function [I, gamma, x_0] = fit_lorentzian(x, y, beta0)
%Function to fit a Lorentzian to a curve using the nlinfit function and
%the function lorentzian.m
%
%[I, gamma, x_0] = fit_lorentzian(x, y, beta0)
%
%INPUTS:
% x     - X values.
% y     - Y values.
% beta0 - Initial values guess. Should be 3 coefficients long,
% corresponding to the I, gamma, and x_0:
% beta0(1)      - x_0
% beta0(2)      - gamma
% beta0(3)      - I
%
beta = nlinfit(x, y, @lorentzian, beta0);
I = beta(3);
gamma = beta(2);
x_0 = beta(1);