#!/bin/sh
# To version bump the Julia code so that I don't have to manually do it myself

cp -r $1 $2
tar -cvf "$1".tar $1/

