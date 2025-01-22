#include "ekf_1.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <eekf/eekf.h>

#include "data_gps.h"
#include "data_imu.h"

#include "usart.h"

#include "MPU9250-DMP.h"
#include "gps_neo6.h"

#define PI 3.14159265
#define G       9.8123

/*==============*/
/* Debug switch */
/*==============*/
#define TEST_DATA

// time step duration
eekf_value dT = 0.01;

eekf_value dT_GPS = 1;
// process noise standard deviation

uint8_t gps_sampling_rate = 1; 
#ifdef TEST_DATA
uint8_t imu_sampling_rate = 100;
#define DURATION_IN_SEC     14
#endif
#ifdef SENSOR_DATA
uint8_t imu_sampling_rate = 50;
#endif

eekf_value std_dev_acc_x = 0.005006432391 * G;
eekf_value std_dev_acc_y = 0.005060815819 * G;
eekf_value std_dev_gyr_z = 0.03773262547;
eekf_value std_dev_gps = 0.01;

uint8_t buffer_ekf_1[128];

uint32_t prediction_time = 0;
uint32_t correction_time = 0;

float a_x_meters = 0.0;
float a_y_meters = 0.0;

void UART2_SendString_ekf_1(char* s)
{
 HAL_UART_Transmit(&huart2, (uint8_t*)s, strlen(s), 1000);
}

eekf_return transition(eekf_mat* xp, eekf_mat* Jf, eekf_mat const *x,
        eekf_mat const *u, void* userData)
{
	EEKF_DECL_MAT_INIT(xu, 4, 1, 0);
	EEKF_DECL_MAT_INIT(B, 4, 2, 0,0,0,0,0,0,0,0);
	
    // the Jacobian of transition() at x
    *EEKF_MAT_EL(*Jf, 0, 0) = 1;
    *EEKF_MAT_EL(*Jf, 1, 0) = 0;
    *EEKF_MAT_EL(*Jf, 2, 0) = 0;
    *EEKF_MAT_EL(*Jf, 3, 0) = 0;
    *EEKF_MAT_EL(*Jf, 0, 1) = 0;
    *EEKF_MAT_EL(*Jf, 1, 1) = 1;
    *EEKF_MAT_EL(*Jf, 2, 1) = 0;
    *EEKF_MAT_EL(*Jf, 3, 1) = 0;
    *EEKF_MAT_EL(*Jf, 0, 2) = dT;
    *EEKF_MAT_EL(*Jf, 1, 2) = 0;
    *EEKF_MAT_EL(*Jf, 2, 2) = 1;
    *EEKF_MAT_EL(*Jf, 3, 2) = 0;
    *EEKF_MAT_EL(*Jf, 0, 3) = 0;
    *EEKF_MAT_EL(*Jf, 1, 3) = dT;
    *EEKF_MAT_EL(*Jf, 2, 3) = 0;
    *EEKF_MAT_EL(*Jf, 3, 3) = 1;

	*EEKF_MAT_EL(B, 0, 0) = dT * dT / 2;
	*EEKF_MAT_EL(B, 1, 0) = 0;
    *EEKF_MAT_EL(B, 2, 0) = dT;
	*EEKF_MAT_EL(B, 3, 0) = 0;
    *EEKF_MAT_EL(B, 0, 1) = 0;
	*EEKF_MAT_EL(B, 1, 1) = dT * dT / 2;
    *EEKF_MAT_EL(B, 2, 1) = 0;
	*EEKF_MAT_EL(B, 3, 1) = dT;
	
    // predict state from current state
    if (NULL == eekf_mat_add(xp, eekf_mat_mul(xp, Jf, x), eekf_mat_mul(&xu, &B, u)))
    {
		printf("arg\n");
        return eEekfReturnComputationFailed;
    }
		
    return eEekfReturnOk;
}

eekf_return measurement(eekf_mat* zp, eekf_mat* Jh, eekf_mat const *x,
        void* userData)
{
    // the Jacobian of measurement() at x
    *EEKF_MAT_EL(*Jh, 0, 0) = 1;
    *EEKF_MAT_EL(*Jh, 1, 0) = 0;
    *EEKF_MAT_EL(*Jh, 0, 1) = 0;
    *EEKF_MAT_EL(*Jh, 1, 1) = 1;

    // compute the measurement from state x
    *EEKF_MAT_EL(*zp, 0, 0) = *EEKF_MAT_EL(*x, 0, 0);
    *EEKF_MAT_EL(*zp, 1, 0) = *EEKF_MAT_EL(*x, 1, 0);

    return eEekfReturnOk;
}

void run_ekf_1(void){
    // filter context
    eekf_context ctx;

    // state of the filter
    #ifdef TEST_DATA
    EEKF_DECL_MAT_INIT(x, 4, 1, data_gps[0][0], data_gps[0][1], 0, 0);
    eekf_value phi = 278.55380 * PI / 180.0;
    eekf_value acc_N = 0.0;
    eekf_value acc_E = 0.0;
    #endif
    EEKF_DECL_MAT_INIT(P, 4, 4, 
		pow(std_dev_acc_x, 2), 0, 0, 0,
		0, pow(std_dev_acc_x, 2), 0, 0,
        0, 0, pow(std_dev_acc_x, 2), 0,
        0, 0, 0, pow(std_dev_acc_x, 2));


    // input and process noise variables            
    EEKF_DECL_MAT_INIT(u, 2, 1, 0, 0);

/*
    EEKF_DECL_MAT_INIT(Q, 4, 4, 
		pow(s_w, 2) * pow(dT, 4) / 4, 0, pow(s_w, 2) * pow(dT, 3) / 2, 0,
		0,  pow(s_w, 2) * pow(dT, 4) / 4, 0, pow(s_w, 2) * pow(dT, 3) / 2,
        pow(s_w, 2) * pow(dT, 3) / 2, 0, pow(s_w, 2) * pow(dT, 2), 0,
        0, 0, 0, pow(s_w, 2) * pow(dT, 2));
*/

    EEKF_DECL_MAT_INIT(Q, 4, 4, 
		pow(std_dev_acc_x, 2) * pow(dT, 2), 0, 0, 0,
		0, pow(std_dev_acc_y, 2) * pow(dT, 2), 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);
        
    // measurement and measurement noise variables
    EEKF_DECL_MAT_INIT(z, 2, 1, 0, 0);
    EEKF_DECL_MAT_INIT(R, 2, 2, pow(std_dev_gps, 2) * pow(dT_GPS, 2), 0,
                                0, pow(std_dev_gps, 2) * pow(dT_GPS, 2));

    eekf_value x_predicted = 0; 
    eekf_value y_predicted = 0; 
    
    #ifdef SENSOR_DATA
    MPU9250_begin();
    MPU9250_setSensors(INV_XYZ_GYRO | INV_XYZ_ACCEL);
    HAL_Delay(500);

    NEO6_Init(&GpsState, &huart1);
    HAL_Delay(500);
    uint32_t Timer_ekf_1 = HAL_GetTick();
    #endif
    
    #ifdef SENSOR_DATA
    /*=====================*/
    /* Initialization data */
    /*=====================*/
    eekf_value phi = 330.0;  /* set at the begining */
    eekf_value acc_x = 0.0;
    eekf_value acc_y = 0.0;
    eekf_value gyr_z = 0.0;
    eekf_value x_0 = 0.0;    /* get data from the GPS */
    eekf_value y_0 = 0.0;    /* get data from the GPS */

    eekf_value acc_x_bias = 0.1049847826;  /* g unit*/
    eekf_value acc_y_bias = 0.02335797101; /* g unit*/
    eekf_value gyr_z_bias = -0.8611884058; 

    eekf_value acc_N = 0.0;
    eekf_value acc_E = 0.0;
    eekf_value x_GPS = 0.0;
    eekf_value y_GPS = 0.0;

    while(1){
        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_ekf_1) > 1000){
            if(NEO6_IsFix(&GpsState)){
                x_0 = GpsState.Latitude;
                y_0 = GpsState.Longitude;
                Timer_ekf_1 = HAL_GetTick();
                break;
            }
        }
        Timer_ekf_1 = HAL_GetTick();
    }
    EEKF_DECL_MAT_INIT(x, 4, 1, x_0, y_0, 0, 0);
    #endif

    // initialize the filter context
    eekf_init(&ctx, &x, &P, transition, measurement, NULL);

    while(1)
    {   
#ifdef TEST_DATA
        for(uint16_t k = 0; k < DURATION_IN_SEC; k++){
            for(uint8_t i = 0; i < imu_sampling_rate/gps_sampling_rate; i++){
                a_x_meters = data_imu[i + k*(imu_sampling_rate/gps_sampling_rate)][4];
                a_y_meters = data_imu[i + k*(imu_sampling_rate/gps_sampling_rate)][5];
                phi += data_imu[i + k*(imu_sampling_rate/gps_sampling_rate)][3];
                if (phi > 2*PI) phi -= 2*PI;
                else if (phi < 0.0) phi += 2*PI;
                acc_N = cos(phi) * a_x_meters - sin(phi) * a_y_meters;
                acc_E = sin(phi) * a_x_meters + cos(phi) * a_y_meters;
                acc_N = acc_N / 111320.0;
                acc_E = acc_E / 111320.0 * cos(*EEKF_MAT_EL(*ctx.x, 0, 0) * PI / 180.0);
                *EEKF_MAT_EL(u, 0, 0) = acc_N;
                *EEKF_MAT_EL(u, 1, 0) = acc_E;
                eekf_predict(&ctx, &u, &Q);   
            }
            *EEKF_MAT_EL(z, 0, 0) = data_gps[k][0];
            *EEKF_MAT_EL(z, 1, 0) = data_gps[k][1];
            x_predicted = *EEKF_MAT_EL(*ctx.x, 0, 0);
            y_predicted = *EEKF_MAT_EL(*ctx.x, 1, 0);
            eekf_correct(&ctx, &z, &R);		
            memset(buffer_ekf_1, 0, 128);
            sprintf((char*)buffer_ekf_1, "N: INS: %.5f   GPS: %.6f   KAL: %.6f   VEL: %.5f\n\r", x_predicted, data_gps[k][0], *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 2, 0));
            UART2_SendString_ekf_1((char*)buffer_ekf_1);
            memset(buffer_ekf_1, 0, 128);
            sprintf((char*)buffer_ekf_1, "E: INS: %.5f   GPS: %.6f   KAL: %.6f   VEL: %.5f\n\r", y_predicted, data_gps[k][1], *EEKF_MAT_EL(*ctx.x, 1, 0), *EEKF_MAT_EL(*ctx.x, 3, 0));
            UART2_SendString_ekf_1((char*)buffer_ekf_1);
        }		
        memset(buffer_ekf_1, 0, 128);
        sprintf((char*)buffer_ekf_1, "==============\n\r");
        UART2_SendString_ekf_1((char*)buffer_ekf_1);
#endif

#ifdef SENSOR_DATA

        if(MPU9250_dataReady()){
            MPU9250_updateAccel();
            MPU9250_updateGyro();
            acc_x = G * MPU9250_calcAccel(ax, acc_x_bias);
            acc_y = G * MPU9250_calcAccel(ay, acc_y_bias);
            gyr_z = MPU9250_calcGyro(gz, gyr_z_bias);
            prediction_time = HAL_GetTick();
            phi += gyr_z;
            if (phi > 360.0) phi -= 360.0;
            else if (phi < 0.0) phi += 360.0;
            phi = (phi / 180.0) * PI;
            acc_N = cos(phi) * acc_x - sin(phi) * acc_y;
            acc_E = sin(phi) * acc_x + cos(phi) * acc_y;
            acc_N = acc_N / 111320.0;
            acc_E = acc_E / 111320.0 * cos(*EEKF_MAT_EL(*ctx.x, 0, 0) * PI / 180.0);
            *EEKF_MAT_EL(u, 0, 0) = acc_N;
            *EEKF_MAT_EL(u, 1, 0) = acc_E;

            eekf_predict(&ctx, &u, &Q);
            prediction_time = HAL_GetTick() - prediction_time;
        } 

        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_ekf_1) > 1000){
            x_predicted = *EEKF_MAT_EL(*ctx.x, 0, 0);
            y_predicted = *EEKF_MAT_EL(*ctx.x, 1, 0);
            if(NEO6_IsFix(&GpsState)){
                x_GPS = GpsState.Latitude;
                y_GPS = GpsState.Longitude; 

                correction_time = HAL_GetTick();
                *EEKF_MAT_EL(z, 0, 0) = x_GPS;
                eekf_correct(&ctx, &z, &R);	
                correction_time = HAL_GetTick() - correction_time;
            }
            else{
                x_GPS = 0.0;
                correction_time = 0.0;
            }
            Timer_ekf_1 = HAL_GetTick();

            memset(buffer_ekf_1, 0, 128);
            sprintf((char*)buffer_ekf_1, "N: INS: %.5f   GPS: %.5f   KAL: %.5f   VEL: %.5f\n\r", x_predicted, x_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 2, 0));
            UART2_SendString_ekf_1((char*)buffer_ekf_1);
            memset(buffer_ekf_1, 0, 128);
            sprintf((char*)buffer_ekf_1, "E: INS: %.5f   GPS: %.5f   KAL: %.5f   VEL: %.5f\n\r", y_predicted, y_GPS, *EEKF_MAT_EL(*ctx.x, 1, 0), *EEKF_MAT_EL(*ctx.x, 3, 0));
            UART2_SendString_ekf_1((char*)buffer_ekf_1);

            /*
            memset(buffer_ekf_1, 0, 128);
            sprintf((char*)buffer_ekf_1, "%.5f, %.5f, %.5f, %.5f\n\r", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0));
            UART2_SendString_ekf_1((char*)buffer_ekf_1);
            */
        }
#endif
    }
}


void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart)
{
	if(huart == GpsState.neo6_huart)
	{
		NEO6_ReceiveUartChar(&GpsState);
	}
}
