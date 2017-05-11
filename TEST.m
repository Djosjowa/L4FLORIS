clear all; close all; clc;
addpath('bin');

a.a = 100;
a.b = 1;


functie = @(x)optimtest(x,a);

options = optimoptions(@fmincon,...
    'Display','iter','Algorithm','interior-point');


[x,fval] = fmincon(functie,[0 0],...
    [],[],[],[],[],[],@unitdisk,options)