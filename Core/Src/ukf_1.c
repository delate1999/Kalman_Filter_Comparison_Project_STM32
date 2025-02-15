#include <string.h>

#include "usart.h"

#include "yafl.h"
#include "ukf_1.h"

#include "data_imu.h"
#include "data_gps.h"

#include "MPU9250-DMP.h"
#include "gps_neo6.h"

/*OK! We can do some filtering...*/

#define PI 3.14159265
#define G       9.8123


#define NX 8 /* State vector dimension       */
#define NZ 3 /* Observation vectio dimention */

#define PI 3.14159265
#define DURATION_IN_SEC     90
#define UART_LOGGING_COUNTER    20

yaflFloat dT_ukf_1 = 0.01;
uint32_t imu_data_index = 0;
uint32_t gps_data_index = 0;
uint32_t uart_index_ukf_1 = 0; 
yaflFloat phi_degree_ukf_1 = 278.55380;
uint8_t buffer_ukf_1[256];

yaflFloat std_dev_acc_x_ukf_1 = 0.005;
yaflFloat std_dev_acc_y_ukf_1 = 0.005;
yaflFloat std_dev_gyr_z_ukf_1 = 0.03;
yaflFloat std_dev_gps_ukf_1 = 0.0005;
yaflFloat std_dev_acc_x_bias_ukf_1 = 0.000001;
yaflFloat std_dev_acc_y_bias_ukf_1 = 0.000001;
yaflFloat std_dev_gyr_z_bias_ukf_1 = 0.000001;

yaflFloat current_acc_x_measured = 0.0;
yaflFloat current_acc_y_measured = 0.0;
yaflFloat current_omega_measured = 0.0;

yaflFloat x_GPS = 0.0;
yaflFloat y_GPS = 0.0;

void UART2_SendString_ukf_1(char* s)
{
 HAL_UART_Transmit(&huart2, (uint8_t*)s, strlen(s), 1000);
}

/* This is our state transition function */
yaflStatusEn fx(yaflKalmanBaseSt * self, yaflFloat * x, yaflFloat * xz)
{
    (void)xz;
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NX == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_2);

    yaflFloat acc_x = current_acc_x_measured - x[5];
    yaflFloat acc_y = current_acc_y_measured - x[6];
    yaflFloat omega = current_omega_measured - x[7];

    if(fabs(omega) <= 0.2){
        omega = 0.0;
    }
    if(fabs(acc_x) <= 0.015 * G){
        acc_x = 0.0;
    }
    if(fabs(acc_y) <= 0.015 * G){
        acc_y = 0.0;
    }

    yaflFloat acc_N = cos((x[2] / 180.0) * PI) * acc_x - sin((x[2] / 180.0) * PI) * acc_y;
    yaflFloat acc_E = sin((x[2] / 180.0) * PI) * acc_x + cos((x[2] / 180.0) * PI) * acc_y;
    acc_N = acc_N / 111320.0;
    acc_E = acc_E / (111320.0 * cos(x[0] * PI / 180.0));

    /*We have a linear uniform motion here:*/
    x[0] = x[0] + x[3] * dT_ukf_1 + (acc_N) * 0.5 * dT_ukf_1 * dT_ukf_1;
    x[1] = x[1] + x[4] * dT_ukf_1 + (acc_E) * 0.5 * dT_ukf_1 * dT_ukf_1;
    x[2] = x[2] + (omega) * dT_ukf_1;
    x[3] = x[3] + (acc_N) * dT_ukf_1;
    x[4] = x[4] + (acc_E) * dT_ukf_1;
    x[5] = x[5];
    x[6] = x[6];
    x[7] = x[7];

    return YAFL_ST_OK;
}

/*This is our measurement function*/
yaflStatusEn hx(yaflKalmanBaseSt * self, yaflFloat * y, yaflFloat * x)
{
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NZ == self->Nz, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(y,             YAFL_ST_INV_ARG_2);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_3);

    /*We see only coordinates, not the velocities*/
    y[0] = x[0];
    y[1] = x[1];
    x[2] = sqrt(x[3] * x[3] + x[4] * x[4]);

    return YAFL_ST_OK;
}

/*---------------------------------------------------------------------------*/
/*This is our filter memory structure*/
typedef struct {
    YAFL_UKF_MEMORY_MIXIN(NX, NZ);
    YAFL_UKF_JULIER_MEMORY_MIXIN(NX, NZ);
    yaflFloat dummy[30];
    /*Other fields*/
} myUKFMemorySt;

/* ==============================*/
/* CALCULATED VALUES AS FOR EKF_2*/
/* ==============================*/
#define DP (0.0000000025)
#define DX_acc (0.0000000025)
#define DX_gyro (0.00000009)
#define DX_bias (1e-14)
#define DZ (0.000025)

myUKFMemorySt ukf_memory =
{
    /*Initial state vector*/
    .x = {
        [0] = 52.3252947,
        [1] = 20.9393912,
        [2] = 245.0,   
        [3] = 0.,
        [4] = 0.,
        [5] = 0.0999847826 * G,
        [6] = 0.02335797101 * G,
        [7] = 1.1061884058
        /*Other values are zeros*/
    },

    /*State covariance components*/
    .Up = {
        /*
        Here we have a unit upper triangular matrix.
        We don't need to store ones, so, only upper parts of three columns are stored
        */
        0,    /*1st column*/
        0,0,  /*2nd column*/
        0,0,0, /*3rd column*/
        0,0,0,0, /*4th column*/
        0,0,0,0,0, /*5th column*/
        0,0,0,0,0,0, /*7th column*/
        0,0,0,0,0,0,0 /*8th column*/
    },
    .Dp = {DP, DP, DP, DP, DP, DP, DP, DP}, /*Diagonal matrix is stored in a vector*/

    /*State noise covariance components*/
    .Uq = {
        0,
        0,0,
        0,0,0,
        0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,0,
        0,0,0,0,0,0,0
    },
    .Dq = {DX_acc, DX_acc, DX_gyro, DX_acc, DX_acc, DX_bias, DX_bias, DX_bias},

    /*Measurement noise covariance components*/
    .Ur = {
        0,
        0, 0
    },
    .Dr = {DZ, DZ, DZ}
};

yaflUKFJulierSt sp = YAFL_UKF_JULIER_INITIALIZER(NX, 0, 0.0, ukf_memory);
yaflUKFSt kf = YAFL_UKF_INITIALIZER(&sp.base, &yafl_ukf_julier_spm, fx, 0, 0, hx, 0, 0, NX, NZ, 0.0, ukf_memory);

/*
Arguments of initializer macro:
fx     - state transition function pointer
jfx    - state transition Jacobian function pointer
hx     - measurement function pointer
jhx    - measurement Jacobian function pointer
zrf    - measurement Residual function pointer (needed to calculate the distance between forecast and measurement vectors)
nx     - the dimension of state vector
nz     - the dimension of measurement vector
rff    - measurement noice covatiance forgetting factor
qff    - process noice covatiance forgetting factor
memory - the name of a memory structure.
*/


void run_ukf_1(void){
    uint32_t prediction_systick_timer = 0; 
    uint32_t correction_systick_timer = 0; 

    yaflStatusEn status;
    yaflFloat z[NZ]; /*Memory for measurement vectors*/

    MPU9250_begin();
    MPU9250_setSensors(INV_XYZ_GYRO | INV_XYZ_ACCEL);
    HAL_Delay(500);

    NEO6_Init(&GpsState, &huart1);
    HAL_Delay(500);
    uint32_t Timer_ukf_1 = HAL_GetTick();

    while(1){
        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_ukf_1) > 1000){
            if(NEO6_IsFix(&GpsState)){
                int x_GPS_int = GpsState.Latitude / 100.0;
                int y_GPS_int = GpsState.Longitude / 100.0;
                x_GPS = x_GPS_int + (GpsState.Latitude - x_GPS_int*100)/60.0;
                y_GPS = y_GPS_int + (GpsState.Longitude - y_GPS_int*100)/60.0; 
                Timer_ukf_1 = HAL_GetTick();
                break;
            }
            Timer_ukf_1 = HAL_GetTick();
        }
    }
    kf.base.base.x[0] = x_GPS;
    kf.base.base.x[1] = y_GPS;

    status = yafl_ukf_post_init(&kf.base);

    while(1){
        if(MPU9250_dataReady()){
            MPU9250_updateAccel();
            MPU9250_updateGyro();
            current_acc_x_measured = G * MPU9250_calcAccel(ax,  0.);
            current_acc_y_measured = (-1.0) * G * MPU9250_calcAccel(ay,  0.);
            current_omega_measured = (-1.0) * MPU9250_calcGyro(gz,  0.);

            
            if(uart_index_ukf_1 >= UART_LOGGING_COUNTER){
                memset(buffer_ukf_1, 0, 256);
                sprintf((char*)buffer_ukf_1, "%.7f, %.7f, %.7f, %.7f, %ld\n", x_GPS, y_GPS, kf.base.base.x[0], kf.base.base.x[1], prediction_systick_timer);
                UART2_SendString_ukf_1((char*)buffer_ukf_1);
                uart_index_ukf_1 = 0; 
            }
            uart_index_ukf_1++;
            
            
            /* This is Bierman filter predict step */
            uint32_t prediction_systick_timer_ms = HAL_GetTick();
            prediction_systick_timer = SysTick->VAL;
            status = yafl_ukf_predict(&kf);
            prediction_systick_timer = prediction_systick_timer - SysTick->VAL + 84000*(HAL_GetTick()-prediction_systick_timer_ms);   
        }                   

        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_ukf_1) > 1000){
            if(NEO6_IsFix(&GpsState)){
                int x_GPS_int = GpsState.Latitude / 100.0;
                int y_GPS_int = GpsState.Longitude / 100.0;
                x_GPS = x_GPS_int + (GpsState.Latitude - x_GPS_int*100)/60.0;
                y_GPS = y_GPS_int + (GpsState.Longitude - y_GPS_int*100)/60.0;

                yaflFloat V_GPS_m_s = fabs(GpsState.SpeedKilometers / 3.6); /* m/s*/
                yaflFloat V_GPS_m_s_n = V_GPS_m_s * cos((kf.base.base.x[2] / 180.0) * PI);
                yaflFloat V_GPS_m_s_e = V_GPS_m_s * sin((kf.base.base.x[2] / 180.0) * PI);
                yaflFloat V_GPS_degree_n = fabs(V_GPS_m_s_n / 111320.0);
                yaflFloat V_GPS_degree_e = fabs(V_GPS_m_s_e / (111320.0 * cos(kf.base.base.x[0] * PI / 180.0)));
                yaflFloat V_GPS_degree_scaled = V_GPS_degree_n * V_GPS_degree_n + V_GPS_degree_e * V_GPS_degree_e;
                yaflFloat V_GPS = sqrt(V_GPS_degree_scaled);

                z[0] = x_GPS;
                z[1] = y_GPS;
                z[2] = V_GPS;

                uint32_t correction_systick_timer_ms = HAL_GetTick();
                correction_systick_timer = SysTick->VAL;
                status = yafl_ukf_update(&kf.base, &z[0]);
                correction_systick_timer = correction_systick_timer - SysTick->VAL + 84000*(HAL_GetTick()-correction_systick_timer_ms);
            }
            else{
               x_GPS = 0.0;
               y_GPS = 0.0;
            }
            Timer_ukf_1 = HAL_GetTick();

            memset(buffer_ukf_1, 0, 256);
            sprintf((char*)buffer_ukf_1, "%.7f, %.7f, %.7f, %.7f, %ld\n", data_gps[gps_data_index][0], data_gps[gps_data_index][1], kf.base.base.x[0], kf.base.base.x[1], correction_systick_timer);
            UART2_SendString_ukf_1((char*)buffer_ukf_1);
        }
        
    }
}