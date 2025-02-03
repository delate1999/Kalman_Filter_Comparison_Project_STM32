#include "ekf_1.h"
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

eekf_value dT_GPS_ekf_1 = 1;

uint8_t gps_sampling_rate_ekf_1 = 1; 
#ifdef TEST_DATA
uint8_t imu_sampling_rate_ekf_1 = 100;
uint8_t imu_data_collected_frequency_ekf_1 = 200;
#define DURATION_IN_SEC     90
#define UART_LOGGING_COUNTER    20
eekf_value dT_ekf_1 = 0.01;
#endif
#ifdef SENSOR_DATA
uint8_t imu_sampling_rate_ekf_1 = 100;
eekf_value dT_ekf_1 = 0.01;
#define UART_LOGGING_COUNTER    20
#endif

uint8_t buffer_ekf_1[256];

void UART2_SendString_ekf_1(char* s)
{
 HAL_UART_Transmit(&huart2, (uint8_t*)s, strlen(s), 1000);
}

eekf_return transition_ekf_1(eekf_mat* xp, eekf_mat* Jf, eekf_mat const *x,
        eekf_mat const *u, void* userData)
{
	// EEKF_DECL_MAT_INIT(xu, 4, 1, 0);
	// EEKF_DECL_MAT_INIT(B, 4, 2, 0,0,0,0,0,0,0,0);
	
    eekf_value x_N = *EEKF_MAT_EL(*x, 0, 0);
    eekf_value y_N = *EEKF_MAT_EL(*x, 1, 0);
    eekf_value phi_degree = *EEKF_MAT_EL(*x, 2, 0);
    eekf_value v = *EEKF_MAT_EL(*x, 4, 0);

    eekf_value a = *EEKF_MAT_EL(*u, 0, 0);
    eekf_value omega = *EEKF_MAT_EL(*u, 1, 0);

    if (phi_degree > 360.0) phi_degree -= 360.0;
    else if (phi_degree < 0.0) phi_degree += 360.0;

    // the Jacobian of transition() at x
    *EEKF_MAT_EL(*Jf, 0, 0) = 1;
    *EEKF_MAT_EL(*Jf, 1, 0) = 0;
    *EEKF_MAT_EL(*Jf, 2, 0) = 0;
    *EEKF_MAT_EL(*Jf, 3, 0) = 0;
    *EEKF_MAT_EL(*Jf, 4, 0) = 0;
    *EEKF_MAT_EL(*Jf, 5, 0) = 0;

    *EEKF_MAT_EL(*Jf, 0, 1) = 0;
    *EEKF_MAT_EL(*Jf, 1, 1) = 1;
    *EEKF_MAT_EL(*Jf, 2, 1) = 0;
    *EEKF_MAT_EL(*Jf, 3, 1) = 0;
    *EEKF_MAT_EL(*Jf, 4, 1) = 0;
    *EEKF_MAT_EL(*Jf, 5, 1) = 0;

    *EEKF_MAT_EL(*Jf, 0, 2) = (-1.0) * sin((phi_degree / 180.0) * PI) * (v * dT_ekf_1 + 0.5 * dT_ekf_1 * dT_ekf_1 * a);
    *EEKF_MAT_EL(*Jf, 1, 2) = cos((phi_degree / 180.0) * PI) * (v * dT_ekf_1 + 0.5 * dT_ekf_1 * dT_ekf_1 * a);
    *EEKF_MAT_EL(*Jf, 2, 2) = 0;
    *EEKF_MAT_EL(*Jf, 3, 2) = 0;
    *EEKF_MAT_EL(*Jf, 4, 2) = 0;
    *EEKF_MAT_EL(*Jf, 5, 2) = 0;

    *EEKF_MAT_EL(*Jf, 0, 3) = 0;
    *EEKF_MAT_EL(*Jf, 1, 3) = 0;
    *EEKF_MAT_EL(*Jf, 2, 3) = dT_ekf_1;
    *EEKF_MAT_EL(*Jf, 3, 3) = 1;
    *EEKF_MAT_EL(*Jf, 4, 3) = 0;
    *EEKF_MAT_EL(*Jf, 5, 3) = 0;

    *EEKF_MAT_EL(*Jf, 0, 4) = cos((phi_degree / 180.0) * PI) * dT_ekf_1;
    *EEKF_MAT_EL(*Jf, 1, 4) = sin((phi_degree / 180.0) * PI) * dT_ekf_1;
    *EEKF_MAT_EL(*Jf, 2, 4) = 0;
    *EEKF_MAT_EL(*Jf, 3, 4) = 0;
    *EEKF_MAT_EL(*Jf, 4, 4) = 1;
    *EEKF_MAT_EL(*Jf, 5, 4) = 0;

    *EEKF_MAT_EL(*Jf, 0, 5) = 0.5 * dT_ekf_1 * dT_ekf_1 * cos((phi_degree / 180.0) * PI);
    *EEKF_MAT_EL(*Jf, 1, 5) = 0.5 * dT_ekf_1 * dT_ekf_1 * sin((phi_degree / 180.0) * PI);
    *EEKF_MAT_EL(*Jf, 2, 5) = 0;
    *EEKF_MAT_EL(*Jf, 3, 5) = 0;
    *EEKF_MAT_EL(*Jf, 4, 5) = dT_ekf_1;
    *EEKF_MAT_EL(*Jf, 5, 5) = 1;

	// *EEKF_MAT_EL(B, 0, 0) = dT * dT / 2;
	// *EEKF_MAT_EL(B, 1, 0) = 0;
    // *EEKF_MAT_EL(B, 2, 0) = dT;
	// *EEKF_MAT_EL(B, 3, 0) = 0;
    // *EEKF_MAT_EL(B, 0, 1) = 0;
	// *EEKF_MAT_EL(B, 1, 1) = dT * dT / 2;
    // *EEKF_MAT_EL(B, 2, 1) = 0;
	// *EEKF_MAT_EL(B, 3, 1) = dT;
	
    *EEKF_MAT_EL(*xp, 0, 0) = x_N + v * dT_ekf_1 * cos((phi_degree / 180.0) * PI) + 0.5 * dT_ekf_1 * dT_ekf_1 * a * cos((phi_degree / 180.0) * PI);
    *EEKF_MAT_EL(*xp, 1, 0) = y_N + v * dT_ekf_1 * sin((phi_degree / 180.0) * PI) + 0.5 * dT_ekf_1 * dT_ekf_1 * a * sin((phi_degree / 180.0) * PI);
    *EEKF_MAT_EL(*xp, 2, 0) = phi_degree + omega * dT_ekf_1;
    *EEKF_MAT_EL(*xp, 3, 0) = omega;
    *EEKF_MAT_EL(*xp, 4, 0) = v + a * dT_ekf_1;
    *EEKF_MAT_EL(*xp, 5, 0) = a;

    // predict state from current state
    // if (NULL == eekf_mat_add(xp, eekf_mat_mul(xp, Jf, x), eekf_mat_mul(&xu, &B, u)))
    // {
	// 	printf("arg\n");
    //     return eEekfReturnComputationFailed;
    // }
		
    return eEekfReturnOk;
}

eekf_return measurement_ekf_1(eekf_mat* zp, eekf_mat* Jh, eekf_mat const *x,
        void* userData)
{
    /* Ver. 1 */
    // the Jacobian of measurement() at x
    // *EEKF_MAT_EL(*Jh, 0, 0) = 1;
    // *EEKF_MAT_EL(*Jh, 1, 0) = 0;
    // *EEKF_MAT_EL(*Jh, 0, 1) = 0;
    // *EEKF_MAT_EL(*Jh, 1, 1) = 1;

    // the Jacobian of measurement() at x
    *EEKF_MAT_EL(*Jh, 0, 0) = 1;
    *EEKF_MAT_EL(*Jh, 1, 0) = 0;
    *EEKF_MAT_EL(*Jh, 0, 1) = 0;
    *EEKF_MAT_EL(*Jh, 1, 1) = 1;
    *EEKF_MAT_EL(*Jh, 0, 2) = 0;
    *EEKF_MAT_EL(*Jh, 1, 2) = 0;
    *EEKF_MAT_EL(*Jh, 0, 3) = 0;
    *EEKF_MAT_EL(*Jh, 1, 3) = 0;
    *EEKF_MAT_EL(*Jh, 0, 4) = 0;
    *EEKF_MAT_EL(*Jh, 1, 4) = 0;
    *EEKF_MAT_EL(*Jh, 0, 5) = 0;
    *EEKF_MAT_EL(*Jh, 1, 5) = 0;

    // compute the measurement from state x
    *EEKF_MAT_EL(*zp, 0, 0) = *EEKF_MAT_EL(*x, 0, 0);
    *EEKF_MAT_EL(*zp, 1, 0) = *EEKF_MAT_EL(*x, 1, 0);

    return eEekfReturnOk;
}

void run_ekf_1(void){
    // filter context
    eekf_context ctx;
    uint32_t uart_index = 0; 
    eekf_value y_GPS = 0.0;
    eekf_value x_GPS = 0.0;
    eekf_value acc_N = 0.0;
    eekf_value acc_E = 0.0;
    eekf_value acc_x = 0.0;
    eekf_value acc_y = 0.0;

    eekf_value acc_geographic = 0.0; 
    eekf_value x_predicted = 0.0;
    eekf_value y_predicted = 0.0;

    uint32_t prediction_time = 0;
    uint32_t correction_time = 0;

    eekf_value std_dev_acc_x = 0.005;
    eekf_value std_dev_acc_y = 0.005;
    eekf_value std_dev_gyr_z = 0.03;
    eekf_value std_dev_gps = 0.0005;

    // state of the filter
    #ifdef TEST_DATA
    uint32_t prediction_systick_timer = 0; 
    uint32_t correction_systick_timer = 0; 
    eekf_value phi = 278.55380;
    eekf_value omega = 0; 

    EEKF_DECL_MAT_INIT(x, 6, 1, data_gps[0][0], data_gps[0][1], phi, 0, 0, 0);
    x_GPS = data_gps[0][0];
    y_GPS = data_gps[0][1];
    #endif
    EEKF_DECL_MAT_INIT(P, 6, 6, 
		pow(std_dev_acc_x, 2) * pow(dT_ekf_1, 2), 0, 0, 0, 0, 0,
		0, pow(std_dev_acc_y, 2) * pow(dT_ekf_1, 2), 0, 0, 0, 0,
        0, 0, pow(std_dev_acc_x, 2) * pow(dT_ekf_1, 2), 0, 0, 0,
        0, 0, 0, pow(std_dev_acc_y, 2)* pow(dT_ekf_1, 2),  0, 0,
        0, 0, 0, 0, pow(std_dev_acc_y, 2)* pow(dT_ekf_1, 2),  0,
        0, 0, 0, 0, 0, pow(std_dev_acc_y, 2)* pow(dT_ekf_1, 2));


    // input and process noise variables            
    EEKF_DECL_MAT_INIT(u, 2, 1, 0, 0);

/*
    EEKF_DECL_MAT_INIT(Q, 4, 4, 
		pow(s_w, 2) * pow(dT, 4) / 4, 0, pow(s_w, 2) * pow(dT, 3) / 2, 0,
		0,  pow(s_w, 2) * pow(dT, 4) / 4, 0, pow(s_w, 2) * pow(dT, 3) / 2,
        pow(s_w, 2) * pow(dT, 3) / 2, 0, pow(s_w, 2) * pow(dT, 2), 0,
        0, 0, 0, pow(s_w, 2) * pow(dT, 2));
*/

    // EEKF_DECL_MAT_INIT(Q, 4, 4, 
	// 	pow(std_dev_acc_x, 2) * pow(dT, 2), 0, 0, 0,
	// 	0, pow(std_dev_acc_y, 2) * pow(dT, 2), 0, 0,
    //     0, 0, pow(std_dev_acc_x, 2) * pow(dT, 2), 0,
    //     0, 0, 0, pow(std_dev_acc_y, 2) * pow(dT, 2));

    EEKF_DECL_MAT_INIT(Q, 6, 6, 
		pow(std_dev_acc_x, 2) * pow(dT_ekf_1, 2), 0, 0, 0, 0, 0,
		0, pow(std_dev_acc_x, 2) * pow(dT_ekf_1, 2), 0, 0, 0, 0,
        0, 0, pow(std_dev_acc_x, 2) * pow(dT_ekf_1, 2), 0, 0, 0,
        0, 0, 0, pow(std_dev_gyr_z, 2) * pow(dT_ekf_1, 2), 0, 0,
        0, 0, 0, 0, pow(std_dev_acc_x, 2) * pow(dT_ekf_1, 2), 0,
        0, 0, 0, 0, 0, pow(std_dev_acc_x, 2) * pow(dT_ekf_1, 2));
        
    // measurement and measurement noise variables
    EEKF_DECL_MAT_INIT(z, 2, 1, 0, 0);
    EEKF_DECL_MAT_INIT(R, 2, 2, pow(std_dev_gps, 2) * pow(dT_GPS_ekf_1, 2), 0,
                                0, pow(std_dev_gps, 2) * pow(dT_GPS_ekf_1, 2));
    
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
    /* set at the begining */
    eekf_value phi_degree = 245.0;   /* Biedronka, Nowodwory */
    eekf_value gyr_z = 0.0;

    eekf_value acc_x_bias = 0.0999847826;  /* g unit*/
    eekf_value acc_y_bias = -0.01335797101; /* g unit*/
    eekf_value gyr_z_bias = -1.1161884058; 

    while(1){
        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_ekf_1) > 1000){
            if(NEO6_IsFix(&GpsState)){
                int x_GPS_int = GpsState.Latitude / 100.0;
                int y_GPS_int = GpsState.Longitude / 100.0;
                x_GPS = x_GPS_int + (GpsState.Latitude - x_GPS_int*100)/60.0;
                y_GPS = y_GPS_int + (GpsState.Longitude - y_GPS_int*100)/60.0; 
                Timer_ekf_1 = HAL_GetTick();
                break;
            }
            Timer_ekf_1 = HAL_GetTick();
        }
    }
    EEKF_DECL_MAT_INIT(x, 6, 1, x_GPS, y_GPS, phi_degree, 0, 0, 0);
    #endif

    eekf_init(&ctx, &x, &P, transition_ekf_1, measurement_ekf_1, NULL);

    while(1)
    {   
#ifdef TEST_DATA
        for(uint16_t k = 0; k < DURATION_IN_SEC; k++){
            for(uint8_t i = 0; i < imu_sampling_rate_ekf_1; i++){
                acc_x = data_imu[(imu_data_collected_frequency_ekf_1/imu_sampling_rate_ekf_1)*i + k*(imu_data_collected_frequency_ekf_1)][1] / dT_ekf_1;
                //acc_y = data_imu[i + k*(imu_sampling_rate)][2];
                //phi += data_imu[i + k*(imu_sampling_rate)][0] * (180.0/PI) * dT;
                omega = data_imu[(imu_data_collected_frequency_ekf_1/imu_sampling_rate_ekf_1)*i + k*(imu_data_collected_frequency_ekf_1)][0] * (180.0/PI) / dT_ekf_1;
                //if (phi > 360.0) phi -= 360.0;
                //else if (phi < 0.0) phi += 360.0;
                phi = *EEKF_MAT_EL(*ctx.x, 2, 0);
                acc_N = cos((phi / 180.0) * PI) * acc_x;// - sin((phi / 180.0) * PI) * acc_y;
                acc_E = sin((phi / 180.0) * PI) * acc_x;// + cos((phi / 180.0) * PI) * acc_y;
                acc_N = acc_N / 111320.0;
                acc_E = acc_E / (111320.0 * cos(*EEKF_MAT_EL(*ctx.x, 0, 0) * PI / 180.0));
                acc_geographic = sqrt(acc_N * acc_N + acc_E * acc_E); 
                //*EEKF_MAT_EL(u, 0, 0) = acc_N;
                //*EEKF_MAT_EL(u, 1, 0) = acc_E;
                *EEKF_MAT_EL(u, 0, 0) = acc_geographic;
                *EEKF_MAT_EL(u, 1, 0) = omega;
                uint32_t prediction_systick_timer_ms = HAL_GetTick();
                prediction_systick_timer = SysTick->VAL;
                eekf_predict(&ctx, &u, &Q);   
                prediction_systick_timer = prediction_systick_timer - SysTick->VAL + 84000*(HAL_GetTick()-prediction_systick_timer_ms);

                
                if(uart_index >= UART_LOGGING_COUNTER){
                    memset(buffer_ekf_1, 0, 256);
                    sprintf((char*)buffer_ekf_1, "%.7f, %.7f, %.7f, %.7f, %ld\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0), prediction_systick_timer);
                    UART2_SendString_ekf_1((char*)buffer_ekf_1);
                    uart_index = 0; 
                }
                uart_index++;
                
                
                
            }
            x_predicted = *EEKF_MAT_EL(*ctx.x, 0, 0);
            y_predicted = *EEKF_MAT_EL(*ctx.x, 1, 0);
            *EEKF_MAT_EL(z, 0, 0) = data_gps[k][0];
            *EEKF_MAT_EL(z, 1, 0) = data_gps[k][1];
            x_GPS = data_gps[k][0];
            y_GPS = data_gps[k][1];
            uint32_t correction_systick_timer_ms = HAL_GetTick();
            correction_systick_timer = SysTick->VAL;
            //if((k <= 54) || (k >= 59)){
            eekf_correct(&ctx, &z, &R);
            //}			
            correction_systick_timer = correction_systick_timer - SysTick->VAL + 84000*(HAL_GetTick()-correction_systick_timer_ms);
            
            /*
            memset(buffer_ekf_1, 0, 256);
            sprintf((char*)buffer_ekf_1, "%.7f, %.7f, %.7f, %.7f, %ld\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0), correction_systick_timer);
            UART2_SendString_ekf_1((char*)buffer_ekf_1);
            */
            /*            
            memset(buffer_ekf_1, 0, 256);
            sprintf((char*)buffer_ekf_1, "N: INS: %.6f   GPS: %.6f   KAL: %.6f   VEL: %.5f\n\r", x_predicted, data_gps[k][0], *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 4, 0));
            UART2_SendString_ekf_1((char*)buffer_ekf_1);
            memset(buffer_ekf_1, 0, 256);
            sprintf((char*)buffer_ekf_1, "E: INS: %.6f   GPS: %.6f   KAL: %.6f   HEAD: %.5f\n\r", y_predicted, data_gps[k][1], *EEKF_MAT_EL(*ctx.x, 1, 0), *EEKF_MAT_EL(*ctx.x, 2, 0));
            UART2_SendString_ekf_1((char*)buffer_ekf_1);
            */
        }	
        
        memset(buffer_ekf_1, 0, 256);
        sprintf((char*)buffer_ekf_1, "==============\n\r");
        UART2_SendString_ekf_1((char*)buffer_ekf_1);
#endif

#ifdef SENSOR_DATA

        if(MPU9250_dataReady()){
            MPU9250_updateAccel();
            MPU9250_updateGyro();
            acc_x = G * MPU9250_calcAccel(ax, acc_x_bias);
            acc_y = (-1.0) * G * MPU9250_calcAccel(ay, acc_y_bias);
            gyr_z = (-1.0) * MPU9250_calcGyro(gz, gyr_z_bias);
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
            //phi_degree -= gyr_z * dT;
            //if (phi_degree > 360.0) phi_degree -= 360.0;
            //else if (phi_degree < 0.0) phi_degree += 360.0;
            //phi_degree = (phi_degree / 180.0) * PI;
            phi_degree = *EEKF_MAT_EL(*ctx.x, 2, 0);
            acc_N = cos((phi_degree / 180.0) * PI) * acc_x;// - sin((phi_degree / 180.0) * PI) * acc_y;
            acc_E = sin((phi_degree / 180.0) * PI) * acc_x;// + cos((phi_degree / 180.0) * PI) * acc_y;
            acc_N = acc_N / 111320.0;
            acc_E = acc_E / (111320.0 * cos(*EEKF_MAT_EL(*ctx.x, 0, 0) * PI / 180.0));
            acc_geographic = sqrt(acc_N * acc_N + acc_E * acc_E);
            *EEKF_MAT_EL(u, 0, 0) = acc_geographic;
            *EEKF_MAT_EL(u, 1, 0) = gyr_z;

            eekf_predict(&ctx, &u, &Q);
            prediction_time = HAL_GetTick() - prediction_time;

            /*
            if(uart_index >= UART_LOGGING_COUNTER){
                memset(buffer_ekf_1, 0, 256);
                sprintf((char*)buffer_ekf_1, "%.7f, %.7f, %.7f, %.7f\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0));
                UART2_SendString_ekf_1((char*)buffer_ekf_1);
                uart_index = 0;
            }
            uart_index++;
            */
            
        } 

        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_ekf_1) > 1000){
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
            Timer_ekf_1 = HAL_GetTick();

            /*
            memset(buffer_ekf_1, 0, 256);
            sprintf((char*)buffer_ekf_1, "N: INS: %.6f   GPS: %.6f   KAL: %.6f   VEL_N: %.5f  ACC_X: %.5f   ACC_N: %.5f   HEAD: %.2f\n\r", x_predicted, x_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 2, 0), acc_x, acc_N, phi_degree);
            UART2_SendString_ekf_1((char*)buffer_ekf_1);
            memset(buffer_ekf_1, 0, 256);
            sprintf((char*)buffer_ekf_1, "E: INS: %.6f   GPS: %.6f   KAL: %.6f   VEL_E: %.5f  ACC_Y: %.5f   ACC_E: %.5f\n\r", y_predicted, y_GPS, *EEKF_MAT_EL(*ctx.x, 1, 0), *EEKF_MAT_EL(*ctx.x, 3, 0), acc_y, acc_E);
            UART2_SendString_ekf_1((char*)buffer_ekf_1);
            */
            
        }

        /*
        memset(buffer_ekf_1, 0, 128);
        sprintf((char*)buffer_ekf_1, "%.5f, %.5f, %.5f, %.5f, %ld, %ld\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0), prediction_time, correction_time);
        UART2_SendString_ekf_1((char*)buffer_ekf_1);
        */
#endif
    }
}

/*
void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart)
{
	if(huart == GpsState.neo6_huart)
	{
		NEO6_ReceiveUartChar(&GpsState);
	}
}
*/


