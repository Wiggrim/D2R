# Doppler2Rotation
An algorithm transforming carrier Doppler & IMU measurements to attitude of the vehicle. Corresponding article is under peer review.

## Features
Being able to work with as few as three satellites.

## Prerequisites
We use [Eigen 3.3.7](http://eigen.tuxfamily.org) for matrix manipulation.

## Example
    mkdir build
    cd build
    cmake ..
    make
    ./dop_rot ../DataSets/data.txt 15 ../outputs/output.txt
