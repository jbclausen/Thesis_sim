#!/bin/bash
# Execute this file to recompile locally
c++ -Wall -shared -fPIC -std=c++11 -O3 -fno-math-errno -fno-trapping-math -ffinite-math-only -I/Users/jensbclausen/opt/miniconda3/envs/fenicstest/include -I/Users/jensbclausen/opt/miniconda3/envs/fenicstest/include/eigen3 -I/Users/jensbclausen/opt/miniconda3/envs/fenicstest/.cache/dijitso/include dolfin_expression_81bcc97c4ff978522a1ef057296cc278.cpp -L/Users/jensbclausen/opt/miniconda3/envs/fenicstest/lib -L/Users/jensbclausen/opt/miniconda3/envs/fenicstest/Users/jensbclausen/opt/miniconda3/envs/fenicstest/lib -L/Users/jensbclausen/opt/miniconda3/envs/fenicstest/.cache/dijitso/lib -Wl,-rpath,/Users/jensbclausen/opt/miniconda3/envs/fenicstest/.cache/dijitso/lib -lpmpi -lmpi -lmpicxx -lpetsc -lslepc -lm -ldl -lz -lpthread -lhdf5 -lboost_timer -ldolfin -Wl,-install_name,/Users/jensbclausen/opt/miniconda3/envs/fenicstest/.cache/dijitso/lib/libdijitso-dolfin_expression_81bcc97c4ff978522a1ef057296cc278.so -olibdijitso-dolfin_expression_81bcc97c4ff978522a1ef057296cc278.so