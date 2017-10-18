#!/bin/bash

for d in */
do
cd $d
pwd
    ${@:1}
cd ../
done
