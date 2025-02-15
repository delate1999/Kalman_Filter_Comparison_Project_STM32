#include "kf_1.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <eekf/eekf.h>

#include "data_gps.h"
#include "data_imu.h"
#include <math.h>

#include "usart.h"

#include "MPU9250-DMP.h"
#include "gps_neo6.h"

#define PI 3.14159265
#define G       9.8123

/*==============*/
/* Debug switch */
/*==============*/
#define TEST_DATA

eekf_value dT_GPS_kf_1 = 1;

uint8_t gps_sampling_rate_kf_1 = 1; 
#ifdef TEST_DATA
uint8_t imu_sampling_rate_kf_1 = 200;
uint8_t imu_data_collected_frequency = 200;
#define DURATION_IN_SEC     90
#define UART_LOGGING_COUNTER    40
eekf_value dT_kf_1 = 0.005;
#endif
#ifdef SENSOR_DATA
uint8_t imu_sampling_rate_kf_1 = 100;
eekf_value dT_kf_1 = 0.01;
#define UART_LOGGING_COUNTER    20
#endif

uint8_t buffer_kf_1[256];

void UART2_SendString_kf_1(char* s)
{
 HAL_UART_Transmit(&huart2, (uint8_t*)s, strlen(s), 1000);
}

eekf_return transition_kf_1(eekf_mat* xp, eekf_mat* Jf, eekf_mat const *x,
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
    *EEKF_MAT_EL(*Jf, 0, 2) = dT_kf_1;
    *EEKF_MAT_EL(*Jf, 1, 2) = 0;
    *EEKF_MAT_EL(*Jf, 2, 2) = 1;
    *EEKF_MAT_EL(*Jf, 3, 2) = 0;
    *EEKF_MAT_EL(*Jf, 0, 3) = 0;
    *EEKF_MAT_EL(*Jf, 1, 3) = dT_kf_1;
    *EEKF_MAT_EL(*Jf, 2, 3) = 0;
    *EEKF_MAT_EL(*Jf, 3, 3) = 1;

    #ifdef TEST_DATA
	*EEKF_MAT_EL(B, 0, 0) = dT_kf_1;
	*EEKF_MAT_EL(B, 1, 0) = 0;
    *EEKF_MAT_EL(B, 2, 0) = 1;
	*EEKF_MAT_EL(B, 3, 0) = 0;
    *EEKF_MAT_EL(B, 0, 1) = 0;
	*EEKF_MAT_EL(B, 1, 1) = dT_kf_1;
    *EEKF_MAT_EL(B, 2, 1) = 0;
	*EEKF_MAT_EL(B, 3, 1) = 1;
    #endif

    #ifdef SENSOR_DATA
	*EEKF_MAT_EL(B, 0, 0) = dT_kf_1 * dT_kf_1 / 2;
	*EEKF_MAT_EL(B, 1, 0) = 0;
    *EEKF_MAT_EL(B, 2, 0) = dT_kf_1;
	*EEKF_MAT_EL(B, 3, 0) = 0;
    *EEKF_MAT_EL(B, 0, 1) = 0;
	*EEKF_MAT_EL(B, 1, 1) = dT_kf_1 * dT_kf_1 / 2;
    *EEKF_MAT_EL(B, 2, 1) = 0;
	*EEKF_MAT_EL(B, 3, 1) = dT_kf_1;
    #endif
	
    // predict state from current state
    if (NULL == eekf_mat_add(xp, eekf_mat_mul(xp, Jf, x), eekf_mat_mul(&xu, &B, u)))
    {
		printf("arg\n");
        return eEekfReturnComputationFailed;
    }
		
    return eEekfReturnOk;
}

eekf_return measurement_kf_1(eekf_mat* zp, eekf_mat* Jh, eekf_mat const *x,
        void* userData)
{
    /* Ver. 1 */
    // the Jacobian of measurement() at x
    // *EEKF_MAT_EL(*Jh, 0, 0) = 1;
    // *EEKF_MAT_EL(*Jh, 1, 0) = 0;
    // *EEKF_MAT_EL(*Jh, 0, 1) = 0;
    // *EEKF_MAT_EL(*Jh, 1, 1) = 1;

    /* Ver. 2*/
    // the Jacobian of measurement() at x
    *EEKF_MAT_EL(*Jh, 0, 0) = 1;
    *EEKF_MAT_EL(*Jh, 1, 0) = 0;

    *EEKF_MAT_EL(*Jh, 0, 1) = 0;
    *EEKF_MAT_EL(*Jh, 1, 1) = 1;

    *EEKF_MAT_EL(*Jh, 0, 2) = 0;
    *EEKF_MAT_EL(*Jh, 1, 2) = 0;

    *EEKF_MAT_EL(*Jh, 0, 3) = 0;
    *EEKF_MAT_EL(*Jh, 1, 3) = 0;

    // compute the measurement from state x
    *EEKF_MAT_EL(*zp, 0, 0) = *EEKF_MAT_EL(*x, 0, 0);
    *EEKF_MAT_EL(*zp, 1, 0) = *EEKF_MAT_EL(*x, 1, 0);

    return eEekfReturnOk;
}

void run_kf_1(void){
    // filter context
    eekf_context ctx;
    uint32_t uart_index = 0; 
    eekf_value y_GPS = 0.0;
    eekf_value x_GPS = 0.0;
    eekf_value acc_N = 0.0;
    eekf_value acc_E = 0.0;
    eekf_value acc_x = 0.0;
    eekf_value acc_y = 0.0;

    uint32_t prediction_time = 0;
    uint32_t correction_time = 0;

    eekf_value std_dev_acc_x = 0.005;
    eekf_value std_dev_acc_y = 0.005;
    eekf_value std_dev_gyr_z = 0.03;
    eekf_value std_dev_gps = 0.001;

    // state of the filter
    #ifdef TEST_DATA
    uint32_t prediction_systick_timer = 0; 
    uint32_t correction_systick_timer = 0; 

    EEKF_DECL_MAT_INIT(x, 4, 1, data_gps[0][0], data_gps[0][1], 0, 0);
    eekf_value phi = 278.55380;
    x_GPS = data_gps[0][0];
    y_GPS = data_gps[0][1];
    #endif
    EEKF_DECL_MAT_INIT(P, 4, 4, 
		pow(std_dev_acc_x, 2) * pow(dT_kf_1, 2), 0, 0, 0,
		0, pow(std_dev_acc_y, 2) * pow(dT_kf_1, 2), 0, 0,
        0, 0, pow(std_dev_acc_x, 2) * pow(dT_kf_1, 2), 0,
        0, 0, 0, pow(std_dev_acc_y, 2)* pow(dT_kf_1, 2));


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
		pow(std_dev_acc_x, 2) * pow(dT_kf_1, 2), 0, 0, 0,
		0, pow(std_dev_acc_y, 2) * pow(dT_kf_1, 2), 0, 0,
        0, 0, pow(std_dev_acc_x, 2) * pow(dT_kf_1, 2), 0,
        0, 0, 0, pow(std_dev_acc_y, 2) * pow(dT_kf_1, 2));
        
    // measurement and measurement noise variables
    EEKF_DECL_MAT_INIT(z, 2, 1, 0, 0);
    EEKF_DECL_MAT_INIT(R, 2, 2, pow(std_dev_gps, 2) * pow(dT_GPS_kf_1, 2), 0,
                                0, pow(std_dev_gps, 2) * pow(dT_GPS_kf_1, 2));
    
    #ifdef SENSOR_DATA
    MPU9250_begin();
    MPU9250_setSensors(INV_XYZ_GYRO | INV_XYZ_ACCEL);
    HAL_Delay(500);

    NEO6_Init(&GpsState, &huart1);
    HAL_Delay(500);
    uint32_t Timer_kf_1 = HAL_GetTick();
    #endif
    
    #ifdef SENSOR_DATA
    /*=====================*/
    /* Initialization data */
    /*=====================*/
    /* set at the begining */
    eekf_value phi_degree = 245.0;   /* Biedronka, Nowodwory */
    eekf_value gyr_z = 0.0;

    eekf_value acc_x_bias = 0.0999847826;  /* g unit*/
    eekf_value acc_y_bias = -0.01335797101; /* g unit*/
    eekf_value gyr_z_bias = -1.0781884058; 

    while(1){
        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_kf_1) > 1000){
            if(NEO6_IsFix(&GpsState)){
                int x_GPS_int = GpsState.Latitude / 100.0;
                int y_GPS_int = GpsState.Longitude / 100.0;
                x_GPS = x_GPS_int + (GpsState.Latitude - x_GPS_int*100)/60.0;
                y_GPS = y_GPS_int + (GpsState.Longitude - y_GPS_int*100)/60.0; 
                Timer_kf_1 = HAL_GetTick();
                break;
            }
            Timer_kf_1 = HAL_GetTick();
        }
    }
    EEKF_DECL_MAT_INIT(x, 4, 1, x_GPS, y_GPS, 0, 0);
    #endif

    eekf_init(&ctx, &x, &P, transition_kf_1, measurement_kf_1, NULL);

    while(1)
    {   
#ifdef TEST_DATA
        for(uint16_t k = 0; k < DURATION_IN_SEC; k++){
            for(uint8_t i = 0; i < imu_sampling_rate_kf_1; i++){
                acc_x = data_imu[(imu_data_collected_frequency/imu_sampling_rate_kf_1)*i + k*(imu_data_collected_frequency)][1];
                acc_y = data_imu[(imu_data_collected_frequency/imu_sampling_rate_kf_1)*i + k*(imu_data_collected_frequency)][2];
                phi += data_imu[(imu_data_collected_frequency/imu_sampling_rate_kf_1)*i + k*(imu_data_collected_frequency)][0] * (180.0/PI);
                if (phi > 360.0) phi -= 360.0;
                else if (phi < 0.0) phi += 360.0;
                acc_N = cos((phi / 180.0) * PI) * acc_x - sin((phi / 180.0) * PI) * acc_y;
                acc_E = sin((phi / 180.0) * PI) * acc_x + cos((phi / 180.0) * PI) * acc_y;
                acc_N = acc_N / 111320.0;
                acc_E = acc_E / (111320.0 * cos(*EEKF_MAT_EL(*ctx.x, 0, 0) * PI / 180.0));
                *EEKF_MAT_EL(u, 0, 0) = acc_N;
                *EEKF_MAT_EL(u, 1, 0) = acc_E;
                prediction_systick_timer = SysTick->VAL;
                eekf_predict(&ctx, &u, &Q);   
                if(SysTick->VAL < prediction_systick_timer)
                    prediction_systick_timer = prediction_systick_timer - SysTick->VAL;
                else 
                    prediction_systick_timer = 0;

                
                if(uart_index >= UART_LOGGING_COUNTER){
                    memset(buffer_kf_1, 0, 256);
                    sprintf((char*)buffer_kf_1, "%.7f, %.7f, %.7f, %.7f, %ld\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0), prediction_systick_timer);
                    UART2_SendString_kf_1((char*)buffer_kf_1);
                    uart_index = 0; 
                }
                uart_index++;
                
                
                
            }
            *EEKF_MAT_EL(z, 0, 0) = data_gps[k][0];
            *EEKF_MAT_EL(z, 1, 0) = data_gps[k][1];
            x_GPS = data_gps[k][0];
            y_GPS = data_gps[k][1];
            correction_systick_timer = SysTick->VAL;
            //if((k <= 54) || (k >= 59)){
            eekf_correct(&ctx, &z, &R);
            //}		
            if(SysTick->VAL < correction_systick_timer)
                correction_systick_timer = correction_systick_timer - SysTick->VAL;
            else 
                correction_systick_timer = 0;
            
            /*
            memset(buffer_kf_1, 0, 256);
            sprintf((char*)buffer_kf_1, "%.7f, %.7f, %.7f, %.7f, %ld\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0), correction_systick_timer);
            UART2_SendString_kf_1((char*)buffer_kf_1);
            */

            /*
            memset(buffer_kf_1, 0, 256);
            sprintf((char*)buffer_kf_1, "N: INS: %.6f   GPS: %.6f   KAL: %.6f   VEL: %.5f\n\r", x_predicted, data_gps[k][0], *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 2, 0));
            UART2_SendString_kf_1((char*)buffer_kf_1);
            memset(buffer_kf_1, 0, 256);
            sprintf((char*)buffer_kf_1, "E: INS: %.6f   GPS: %.6f   KAL: %.6f   VEL: %.5f\n\r", y_predicted, data_gps[k][1], *EEKF_MAT_EL(*ctx.x, 1, 0), *EEKF_MAT_EL(*ctx.x, 3, 0));
            UART2_SendString_kf_1((char*)buffer_kf_1);
            */
        }	
        
        memset(buffer_kf_1, 0, 256);
        sprintf((char*)buffer_kf_1, "==============\n\r");
        UART2_SendString_kf_1((char*)buffer_kf_1);
#endif

#ifdef SENSOR_DATA

        if(MPU9250_dataReady()){
            MPU9250_updateAccel();
            MPU9250_updateGyro();
            acc_x = G * MPU9250_calcAccel(ax, acc_x_bias);
            acc_y = (-1.0) * G * MPU9250_calcAccel(ay, acc_y_bias);
            gyr_z = MPU9250_calcGyro(gz, gyr_z_bias);
            prediction_time = HAL_GetTick();
            if(fabs(gyr_z) <= 0.2){
                gyr_z = 0.0;
            }
            if(fabs(acc_x) <= 0.015 * G){
                acc_x = 0.0;
            }
            if(fabs(acc_y) <= 0.015 * G){
                acc_y = 0.0;
            }
            phi_degree -= gyr_z * dT_kf_1;
            if (phi_degree > 360.0) phi_degree -= 360.0;
            else if (phi_degree < 0.0) phi_degree += 360.0;
            //phi_degree = (phi_degree / 180.0) * PI;
            acc_N = cos((phi_degree / 180.0) * PI) * acc_x - sin((phi_degree / 180.0) * PI) * acc_y;
            acc_E = sin((phi_degree / 180.0) * PI) * acc_x + cos((phi_degree / 180.0) * PI) * acc_y;
            acc_N = acc_N / 111320.0;
            acc_E = acc_E / (111320.0 * cos(*EEKF_MAT_EL(*ctx.x, 0, 0) * PI / 180.0));
            *EEKF_MAT_EL(u, 0, 0) = acc_N;
            *EEKF_MAT_EL(u, 1, 0) = acc_E;

            eekf_predict(&ctx, &u, &Q);
            prediction_time = HAL_GetTick() - prediction_time;

            /*
            if(uart_index >= UART_LOGGING_COUNTER){
                memset(buffer_kf_1, 0, 256);
                sprintf((char*)buffer_kf_1, "%.7f, %.7f, %.7f, %.7f\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0));
                UART2_SendString_kf_1((char*)buffer_kf_1);
                uart_index = 0;
            }
            uart_index++;
            */
        } 

        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_kf_1) > 1000){
            if(NEO6_IsFix(&GpsState)){
                int x_GPS_int = GpsState.Latitude / 100.0;
                int y_GPS_int = GpsState.Longitude / 100.0;
                x_GPS = x_GPS_int + (GpsState.Latitude - x_GPS_int*100)/60.0;
                y_GPS = y_GPS_int + (GpsState.Longitude - y_GPS_int*100)/60.0;     

                correction_time = HAL_GetTick();
                *EEKF_MAT_EL(z, 0, 0) = x_GPS;
                *EEKF_MAT_EL(z, 1, 0) = y_GPS;
                eekf_correct(&ctx, &z, &R);	
                correction_time = HAL_GetTick() - correction_time;
            }
            else{
                x_GPS = 0.0;
                y_GPS = 0.0;
                correction_time = 0.0;
            }
            Timer_kf_1 = HAL_GetTick();

            /*
            memset(buffer_kf_1, 0, 256);
            sprintf((char*)buffer_kf_1, "N: INS: %.6f   GPS: %.6f   KAL: %.6f   VEL_N: %.5f  ACC_X: %.5f   ACC_N: %.5f   HEAD: %.2f\n\r", x_predicted, x_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 2, 0), acc_x, acc_N, phi_degree);
            UART2_SendString_kf_1((char*)buffer_kf_1);
            memset(buffer_kf_1, 0, 256);
            sprintf((char*)buffer_kf_1, "E: INS: %.6f   GPS: %.6f   KAL: %.6f   VEL_E: %.5f  ACC_Y: %.5f   ACC_E: %.5f\n\r", y_predicted, y_GPS, *EEKF_MAT_EL(*ctx.x, 1, 0), *EEKF_MAT_EL(*ctx.x, 3, 0), acc_y, acc_E);
            UART2_SendString_kf_1((char*)buffer_kf_1);
            */
            
        }

        
        memset(buffer_kf_1, 0, 128);
        sprintf((char*)buffer_kf_1, "%.5f, %.5f, %.5f, %.5f, %d, %d\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0), prediction_time, correction_time);
        UART2_SendString_kf_1((char*)buffer_kf_1);
    
#endif
    }
}


