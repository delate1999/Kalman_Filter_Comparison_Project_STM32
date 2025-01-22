/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2025 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */ 
 #define EKF_1
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "i2c.h"
#include "usart.h"
#include "gpio.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#ifdef IMU
#include "mpu6050.h"
#include "MPU9250-DMP.h"
#endif

#ifdef GPS
#include "gps_neo6.h"
#endif

#ifdef EKF_1
#include "ekf_1.h"
#endif

#include "data_gps.h"
#include "data_imu.h"
#include <string.h>
#include <stdio.h>
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */
#ifdef TEST
uint8_t buffer[128];
#endif

#ifdef IMU 
float acc_x, acc_y, acc_z, gyr_x, gyr_y, gyr_z;
float mag_x, mag_y, mag_z;
uint8_t buffer[128];
float magnetic_declination = 6.83;
float acc_x_bias = 0.1049847826;  /* g unit*/
float acc_y_bias = 0.02335797101; /* g unit*/
float gyr_z_bias = -0.8611884058; /* g unit*/
#endif

#ifdef GPS
uint8_t Message[64];
uint8_t MessageLength;
#endif

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

void UART2_SendString(char* s)
{
 HAL_UART_Transmit(&huart2, (uint8_t*)s, strlen(s), 1000);
}

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_USART2_UART_Init();
  MX_I2C1_Init();
  MX_USART1_UART_Init();
  /* USER CODE BEGIN 2 */

  #ifdef IMU 
  MPU9250_begin();
  MPU9250_setSensors(INV_XYZ_GYRO | INV_XYZ_ACCEL | INV_XYZ_COMPASS);
  //MPU6050_Init(&hi2c1);
  HAL_Delay(500);
  #endif

  #ifdef GPS
  NEO6_Init(&GpsState, &huart1);
  HAL_Delay(500);
  uint32_t Timer = HAL_GetTick();
  #endif

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  { 
    #ifdef IMU   
    //MPU6050_GetAccelerometerScaled(&ax, &ay, &az);
	  //MPU6050_GetGyroscopeScaled(&gx, &gy, &gz);
    //AK8963_GetMagnetometerScaled(&mag_x, &mag_y, &mag_z);
    if(MPU9250_dataReady()){
      MPU9250_updateCompass();
      MPU9250_updateAccel();
      MPU9250_updateGyro();
      acc_x = MPU9250_calcAccel(ax, acc_x_bias);
      acc_y = MPU9250_calcAccel(ay, acc_y_bias);
      acc_z = MPU9250_calcAccel(az, 0.0);
      gyr_x = MPU9250_calcGyro(gx, 0.0);
      gyr_y = MPU9250_calcGyro(gy, 0.0);
      gyr_z = MPU9250_calcGyro(gz, gyr_z_bias);
      heading = MPU9250_computeCompassHeading();
      heading += magnetic_declination;
      if (heading > 360.0) heading -= 360.0;
      else if (heading < 0.0) heading += 360.0;
      memset(buffer, 0, 128);
      sprintf((char*)buffer, "%.4f,%.4f,%.4f\n\r", acc_x, acc_y, gyr_z);
      UART2_SendString((char*)buffer);
    }
    HAL_Delay(200);
    #endif

    #ifdef GPS
    NEO6_Task(&GpsState);

	  if((HAL_GetTick() - Timer) > 1000)
	  {
		  MessageLength = sprintf((char*)Message, "\033[2J\033[;H"); // Clear terminal and home cursor
		  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

		  if(NEO6_IsFix(&GpsState))
		  {
			  MessageLength = sprintf((char*)Message, "UTC Time: %02d:%02d:%02d\n\r", GpsState.Hour, GpsState.Minute, GpsState.Second);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

			  MessageLength = sprintf((char*)Message, "Date: %02d.%02d.20%02d\n\r", GpsState.Day, GpsState.Month, GpsState.Year);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

			  MessageLength = sprintf((char*)Message, "Latitude: %.2f %c\n\r", GpsState.Latitude, GpsState.LatitudeDirection);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

			  MessageLength = sprintf((char*)Message, "Longitude: %.2f %c\n\r", GpsState.Longitude, GpsState.LongitudeDirection);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

			  MessageLength = sprintf((char*)Message, "Altitude: %.2f m above sea level\n\r", GpsState.Altitude);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

			  MessageLength = sprintf((char*)Message, "Speed: %.2f knots, %f km/h\n\r", GpsState.SpeedKnots, GpsState.SpeedKilometers);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

			  MessageLength = sprintf((char*)Message, "Satelites: %d\n\r", GpsState.SatelitesNumber);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

			  MessageLength = sprintf((char*)Message, "Dilution of precision: %.2f\n\r", GpsState.Dop);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

			  MessageLength = sprintf((char*)Message, "Horizontal dilution of precision: %.2f\n\r", GpsState.Hdop);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);

			  MessageLength = sprintf((char*)Message, "Vertical dilution of precision: %.2f\n\r", GpsState.Vdop);
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);
		  }
		  else
		  {
			  MessageLength = sprintf((char*)Message, "No Fix\n\r");
			  HAL_UART_Transmit(&huart2, Message, MessageLength, 1000);
		  }

		  Timer = HAL_GetTick();
	  }
    #endif

    #ifdef TEST
    /*
    float delta_t_GPS = 1.0;
    float delta_t_IMU = 0.01;
    float x_old = data_gps[0][0]; 
    float v_old = 0.0; 
    float phi_old = 278.55380 * 3.14159265 / 180.0;
    float Q = 100;
    float P = 0;
    float R = 9; 

    float x_predicted, v_predicted, x_kalman, K;
    uint16_t index = 0; 
    for(uint16_t k = 0; k < 4; k++){
      for(uint16_t i = 0; i < 100; i++){
        x_predicted = x_old +  v_old*delta_t_IMU + 0.5*delta_t_IMU*delta_t_IMU*data_imu[index][3];
        v_predicted = v_old + delta_t_IMU*data_imu[index][3];
        x_old = x_predicted;
        v_old = v_predicted;
        index++;
      }
      P = P + Q;
      K = P*(1/(P + R));
      x_kalman = x_predicted + K * (data_gps[k][0] - x_predicted);
      P = (1 - K)/P;

      memset(buffer, 0, 128);
      sprintf((char*)buffer, "INS: %.4f   GPS: %.4f   KAL: %.4f\n\r", x_predicted, data_gps[k][0], x_kalman);
      UART2_SendString((char*)buffer);
    }
    memset(buffer, 0, 128);
    */
    #endif

    #ifdef EKF_1
    run_ekf_1();
    #endif
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE3);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLM = 16;
  RCC_OscInitStruct.PLL.PLLN = 336;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV4;
  RCC_OscInitStruct.PLL.PLLQ = 2;
  RCC_OscInitStruct.PLL.PLLR = 2;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
  {
    Error_Handler();
  }
}

#if defined GPS
/* USER CODE BEGIN 4 */
void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart)
{
	if(huart == GpsState.neo6_huart)
	{
		NEO6_ReceiveUartChar(&GpsState);
	}
}
#endif
/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
