#include "ekf_2.h"
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

eekf_value dT_GPS = 1;

uint8_t gps_sampling_rate = 1; 
uint8_t imu_sampling_rate = 100;
eekf_value dT = 0.01;
#define UART_LOGGING_COUNTER    20

eekf_value std_dev_acc_x = 0.005;
eekf_value std_dev_acc_y = 0.005;
eekf_value std_dev_gyr_z = 0.03;
eekf_value std_dev_gps = 0.0005;

uint8_t buffer_ekf_2[256];

void UART2_SendString_ekf_2(char* s)
{
 HAL_UART_Transmit(&huart2, (uint8_t*)s, strlen(s), 1000);
}

eekf_return transition(eekf_mat* xp, eekf_mat* Jf, eekf_mat const *x,
        eekf_mat const *u, void* userData)
{
    EEKF_DECL_MAT_INIT(xu, 5, 1, 0);
	EEKF_DECL_MAT_INIT(B, 5, 3, 0);
	
    // the Jacobian of transition() at x
    *EEKF_MAT_EL(*Jf, 0, 0) = 1;
    *EEKF_MAT_EL(*Jf, 1, 0) = 0;
    *EEKF_MAT_EL(*Jf, 2, 0) = 0;
    *EEKF_MAT_EL(*Jf, 3, 0) = 0;
    *EEKF_MAT_EL(*Jf, 4, 0) = 0;

    *EEKF_MAT_EL(*Jf, 0, 1) = 0;
    *EEKF_MAT_EL(*Jf, 1, 1) = 1;
    *EEKF_MAT_EL(*Jf, 2, 1) = 0;
    *EEKF_MAT_EL(*Jf, 3, 1) = 0;
    *EEKF_MAT_EL(*Jf, 4, 1) = 0;

    *EEKF_MAT_EL(*Jf, 0, 2) = 0;
    *EEKF_MAT_EL(*Jf, 1, 2) = 0;
    *EEKF_MAT_EL(*Jf, 2, 2) = 1;
    *EEKF_MAT_EL(*Jf, 3, 2) = 0;
    *EEKF_MAT_EL(*Jf, 4, 2) = 0;

    *EEKF_MAT_EL(*Jf, 0, 3) = dT;
    *EEKF_MAT_EL(*Jf, 1, 3) = 0;
    *EEKF_MAT_EL(*Jf, 2, 3) = 0;
    *EEKF_MAT_EL(*Jf, 3, 3) = 1;
    *EEKF_MAT_EL(*Jf, 4, 3) = 0;

    *EEKF_MAT_EL(*Jf, 0, 4) = 0;
    *EEKF_MAT_EL(*Jf, 1, 4) = dT;
    *EEKF_MAT_EL(*Jf, 2, 4) = 0;
    *EEKF_MAT_EL(*Jf, 3, 4) = 0;
    *EEKF_MAT_EL(*Jf, 4, 4) = 1;


	*EEKF_MAT_EL(B, 0, 0) = dT * dT / 2;
	*EEKF_MAT_EL(B, 1, 0) = 0;
    *EEKF_MAT_EL(B, 2, 0) = 0;
	*EEKF_MAT_EL(B, 3, 0) = dT;
    *EEKF_MAT_EL(B, 4, 0) = 0;

    *EEKF_MAT_EL(B, 0, 1) = 0;
	*EEKF_MAT_EL(B, 1, 1) = dT * dT / 2;
    *EEKF_MAT_EL(B, 2, 1) = 0;
	*EEKF_MAT_EL(B, 3, 1) = 0;
    *EEKF_MAT_EL(B, 4, 1) = dT;

    *EEKF_MAT_EL(B, 0, 2) = 0;
	*EEKF_MAT_EL(B, 1, 2) = 0;
    *EEKF_MAT_EL(B, 2, 2) = dT;
	*EEKF_MAT_EL(B, 3, 2) = 0;
    *EEKF_MAT_EL(B, 4, 2) = 0;
	
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

    eekf_value V_x = *EEKF_MAT_EL(*x, 3, 0);
    eekf_value V_y = *EEKF_MAT_EL(*x, 4, 0);

    // the Jacobian of measurement() at x
    *EEKF_MAT_EL(*Jh, 0, 0) = 1;
    *EEKF_MAT_EL(*Jh, 1, 0) = 0;
    *EEKF_MAT_EL(*Jh, 2, 0) = 0;

    *EEKF_MAT_EL(*Jh, 0, 1) = 0;
    *EEKF_MAT_EL(*Jh, 1, 1) = 1;
    *EEKF_MAT_EL(*Jh, 2, 1) = 0;

    *EEKF_MAT_EL(*Jh, 0, 2) = 0;
    *EEKF_MAT_EL(*Jh, 1, 2) = 0;
    *EEKF_MAT_EL(*Jh, 2, 2) = 0;

    *EEKF_MAT_EL(*Jh, 0, 3) = 0;
    *EEKF_MAT_EL(*Jh, 1, 3) = 0;
    *EEKF_MAT_EL(*Jh, 2, 3) = V_x / sqrt(V_x * V_x + V_y * V_y);

    *EEKF_MAT_EL(*Jh, 0, 4) = 0;
    *EEKF_MAT_EL(*Jh, 1, 4) = 0;
    *EEKF_MAT_EL(*Jh, 2, 4) = V_y / sqrt(V_x * V_x + V_y * V_y);

    // compute the measurement from state x
    *EEKF_MAT_EL(*zp, 0, 0) = *EEKF_MAT_EL(*x, 0, 0);
    *EEKF_MAT_EL(*zp, 1, 0) = *EEKF_MAT_EL(*x, 1, 0);
    *EEKF_MAT_EL(*zp, 2, 0) = sqrt(V_x * V_x + V_y * V_y);

    return eEekfReturnOk;
}

void run_ekf_2(void){

    uint32_t prediction_systick_timer = 0; 
    uint32_t correction_systick_timer = 0; 

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

    EEKF_DECL_MAT_INIT(P, 5, 5, 
		pow(std_dev_acc_x, 2) * pow(dT, 2), 0, 0, 0, 0,
		0, pow(std_dev_acc_y, 2) * pow(dT, 2), 0, 0, 0,
        0, 0, pow(std_dev_gyr_z, 2) * pow(dT, 2), 0, 0,
        0, 0, 0, pow(std_dev_acc_x, 2)* pow(dT, 2),  0,
        0, 0, 0, 0, pow(std_dev_acc_y, 2)* pow(dT, 2));


    // input and process noise variables            
    EEKF_DECL_MAT_INIT(u, 3, 1, 0);

    EEKF_DECL_MAT_INIT(Q, 5, 5, 
		pow(std_dev_acc_x, 2) * pow(dT, 2), 0, 0, 0, 0,
		0, pow(std_dev_acc_y, 2) * pow(dT, 2), 0, 0, 0,
        0, 0, pow(std_dev_gyr_z, 2) * pow(dT, 2), 0, 0,
        0, 0, 0, pow(std_dev_gyr_z, 2) * pow(dT, 2), 0,
        0, 0, 0, 0, pow(std_dev_acc_y, 2) * pow(dT, 2));
        
    // measurement and measurement noise variables
    EEKF_DECL_MAT_INIT(z, 3, 1, 0, 0, 0);
    EEKF_DECL_MAT_INIT(R, 3, 3, pow(std_dev_gps, 2) * pow(dT_GPS, 2), 0, 0,
                                0, pow(std_dev_gps, 2) * pow(dT_GPS, 2), 0,
                                0, 0, pow(std_dev_gps, 2) * pow(dT_GPS, 2));
    
    MPU9250_begin();
    MPU9250_setSensors(INV_XYZ_GYRO | INV_XYZ_ACCEL);
    HAL_Delay(500);

    NEO6_Init(&GpsState, &huart1);
    HAL_Delay(500);
    uint32_t Timer_ekf_2 = HAL_GetTick();
    
    /*=====================*/
    /* Initialization data */
    /*=====================*/
    /* set at the begining */
    eekf_value phi_degree = 245.0;   /* Biedronka, Nowodwory */
    eekf_value gyr_z = 0.0;

    eekf_value acc_x_bias = 0.0999847826;  /* g unit*/
    eekf_value acc_y_bias = -0.01335797101; /* g unit*/
    eekf_value gyr_z_bias = -1.1061884058; 

    while(1){
        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_ekf_2) > 1000){
            if(NEO6_IsFix(&GpsState)){
                int x_GPS_int = GpsState.Latitude / 100.0;
                int y_GPS_int = GpsState.Longitude / 100.0;
                x_GPS = x_GPS_int + (GpsState.Latitude - x_GPS_int*100)/60.0;
                y_GPS = y_GPS_int + (GpsState.Longitude - y_GPS_int*100)/60.0; 
                Timer_ekf_2 = HAL_GetTick();
                break;
            }
            Timer_ekf_2 = HAL_GetTick();
        }
    }

    EEKF_DECL_MAT_INIT(x, 5, 1, x_GPS, y_GPS, phi_degree, 0, 0);

    eekf_init(&ctx, &x, &P, transition, measurement, NULL);

    while(1)
    {   
        if(MPU9250_dataReady()){
            MPU9250_updateAccel();
            MPU9250_updateGyro();
            acc_x = G * MPU9250_calcAccel(ax, acc_x_bias);
            acc_y = (-1.0) * G * MPU9250_calcAccel(ay, acc_y_bias);
            gyr_z = (-1.0) * MPU9250_calcGyro(gz, gyr_z_bias);
            if(fabs(gyr_z) <= 0.2){
                gyr_z = 0.0;
            }
            if(fabs(acc_x) <= 0.02 * G){
                acc_x = 0.0;
            }
            if(fabs(acc_y) <= 0.03 * G){
                acc_y = 0.0;
            }
            phi_degree = *EEKF_MAT_EL(*ctx.x, 2, 0);
            if (phi_degree > 360.0) phi_degree -= 360.0;
            else if (phi_degree < 0.0) phi_degree += 360.0;
            acc_N = cos((phi_degree / 180.0) * PI) * acc_x - sin((phi_degree / 180.0) * PI) * acc_y;
            acc_E = sin((phi_degree / 180.0) * PI) * acc_x + cos((phi_degree / 180.0) * PI) * acc_y;
            acc_N = acc_N / 111320.0;
            acc_E = acc_E / (111320.0 * cos(*EEKF_MAT_EL(*ctx.x, 0, 0) * PI / 180.0));
            *EEKF_MAT_EL(u, 0, 0) = acc_N;
            *EEKF_MAT_EL(u, 1, 0) = acc_E;
            *EEKF_MAT_EL(u, 2, 0) = gyr_z;
            uint32_t prediction_systick_timer_ms = HAL_GetTick();
            prediction_systick_timer = SysTick->VAL;
            eekf_predict(&ctx, &u, &Q);
            prediction_systick_timer = prediction_systick_timer - SysTick->VAL + 84000*(HAL_GetTick()-prediction_systick_timer_ms);

            if(uart_index >= UART_LOGGING_COUNTER){
                memset(buffer_ekf_2, 0, 256);
                sprintf((char*)buffer_ekf_2, "%.7f, %.7f, %.7f, %.7f, %ld\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0), prediction_systick_timer);
                UART2_SendString_ekf_2((char*)buffer_ekf_2);
                uart_index = 0;
            }
            uart_index++;
        } 

        NEO6_Task(&GpsState);
        if((HAL_GetTick() - Timer_ekf_2) > 1000){
            if(NEO6_IsFix(&GpsState)){
                int x_GPS_int = GpsState.Latitude / 100.0;
                int y_GPS_int = GpsState.Longitude / 100.0;
                x_GPS = x_GPS_int + (GpsState.Latitude - x_GPS_int*100)/60.0;
                y_GPS = y_GPS_int + (GpsState.Longitude - y_GPS_int*100)/60.0;

                phi_degree = *EEKF_MAT_EL(*ctx.x, 2, 0);
                eekf_value V_GPS_m_s = fabs(GpsState.SpeedKilometers / 3.6); /* m/s*/
                eekf_value V_GPS_m_s_n = V_GPS_m_s * cos((phi_degree / 180.0) * PI);
                eekf_value V_GPS_m_s_e = V_GPS_m_s * sin((phi_degree / 180.0) * PI);
                eekf_value V_GPS_degree_n = fabs(V_GPS_m_s_n / 111320.0);
                eekf_value V_GPS_degree_e = fabs(V_GPS_m_s_e / (111320.0 * cos(*EEKF_MAT_EL(*ctx.x, 0, 0) * PI / 180.0)));
                eekf_value V_GPS_degree_scaled = V_GPS_degree_n * V_GPS_degree_n + V_GPS_degree_e * V_GPS_degree_e;
                eekf_value V_GPS = sqrt(V_GPS_degree_scaled);

                *EEKF_MAT_EL(z, 0, 0) = x_GPS;
                *EEKF_MAT_EL(z, 1, 0) = y_GPS;
                *EEKF_MAT_EL(z, 2, 0) = V_GPS;
                uint32_t correction_systick_timer_ms = HAL_GetTick();
                correction_systick_timer = SysTick->VAL;
                eekf_correct(&ctx, &z, &R);	
                correction_systick_timer = correction_systick_timer - SysTick->VAL + 84000*(HAL_GetTick()-correction_systick_timer_ms);
            }
            else{
               x_GPS = 0.0;
               y_GPS = 0.0;
            }
            Timer_ekf_2 = HAL_GetTick();
            
            memset(buffer_ekf_2, 0, 128);
            sprintf((char*)buffer_ekf_2, "%.5f, %.5f, %.5f, %.5f, %ld\n", x_GPS, y_GPS, *EEKF_MAT_EL(*ctx.x, 0, 0), *EEKF_MAT_EL(*ctx.x, 1, 0), correction_systick_timer);
            UART2_SendString_ekf_2((char*)buffer_ekf_2);
            
        }
    }
}
