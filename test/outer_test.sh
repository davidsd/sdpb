#!/bin/bash

./build/outer_limits --sdp test/toy_functions.json --points test/toy_functions_points.json --precision=128 --dualityGapThreshold=1e-10 --primalErrorThreshold=1e-10 --dualErrorThreshold=1e-10 --initialMatrixScalePrimal=1e1 --initialMatrixScaleDual=1e1 --maxIterations=1000 --verbosity=0
