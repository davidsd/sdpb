#!/bin/bash

./build/spectrum --precision=128 -i test/spectrum_test.json --solution test/spectrum_test_out/y.txt -o spectrum_test_spectrum.json --threshold=1e-20 --format=PMP
