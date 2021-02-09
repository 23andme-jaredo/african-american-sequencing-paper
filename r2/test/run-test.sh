#! /usr/bin/env bash

bin/r2 guess.bcf truth.bcf

bin/r2 --missing-as-homref guess.bcf truth.bcf ##Note there are no indels on the chip so VAR[GUESS]==0 hence covariance is undefined
