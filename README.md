# Project Title

Kalman Filter implementations comparison on ARM MCU. 

## Description

Project compares KF, EKF and UKF in 2D position estimation in GPS/INS integration. Sensors MPU9250 and NEO-7M were used, converted data from 
awesome-gins-datasets were also included to work without sensors. 

In _papers_ folder there is thesis which used this repository. The paper might better explain usage and concept behind repository.

### Dependencies

* STM32CubeMX
* STM32F44x

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-gins-datasets](https://github.com/i2Nav-WHU/awesome-gins-datasets/tree/main)
* [eekf](https://github.com/dr-duplo/eekf/tree/master)
* [yafl](https://github.com/shkolnick-kun/yafl)