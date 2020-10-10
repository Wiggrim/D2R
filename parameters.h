#pragma once

const double DOP_STD = 5e-2;
const double OSC_STD = 5e-2;

const double INIT_GYR_STD = 0.2 / 180 * 3.14159;
const double GYR_NOISE_STD = 1e-3;
const double GYR_DRIFT_STD = 1e-5;

const double INIT_ACC_STD = 5e-2;
const double ACC_NOISE_STD = 1e-3;
const double ACC_DRIFT_STD = 1e-5;

const double CARRIER_FREQ = 1575.42e6;
const bool INVERSE_DOP = false;

