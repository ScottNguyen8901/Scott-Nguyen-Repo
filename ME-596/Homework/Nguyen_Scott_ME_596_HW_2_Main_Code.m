close all
clear
clc

syms ct1 st1 ct2 st2 ct3 st3

R_1_t1 = [1 0 0;
       0 ct1 st1;
       0 -st1 ct1];

R_2_t2 = [ct2 0 -st2;
          0 1 0;
          st2 0 ct2];

R_1_t3 = [1 0 0;
       0 ct3 st3;
       0 -st3 ct3];

R_121 = R_1_t1 * R_2_t2 * R_1_t3
