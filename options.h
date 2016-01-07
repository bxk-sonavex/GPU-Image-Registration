/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org / sf.chen@ieee.org
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#ifndef FLIRT_OPTIONS_H
#define FLIRT_OPTIONS_H

#define FLIRT_COST_FUNC_MI 0
#define FLIRT_COST_FUNC_CR 1
#define FLIRT_COST_FUNC_NMI 2

//8mm stage
#define FLIRT_COARSE_ANGLE_STEP			+12.0f
#define FLIRT_FINER_ANGLE_STEP			(FLIRT_COARSE_ANGLE_STEP/2.0f)
#define FLIRT_ANGLE_MIN					-36.0f
#define FLIRT_ANGLE_MAX					+36.0f
#define FLIRT_TRANS_XY_STEP				+4.0f
#define FLIRT_TRANS_XY_MIN				-8.0f
#define FLIRT_TRANS_XY_MAX				+8.0f
#define FLIRT_TRANS_Z_STEP				+4.0f
#define FLIRT_TRANS_Z_MIN				-20.0f
#define FLIRT_TRANS_Z_MAX				+20.0f
#define FLIRT_SCALE_STEP				+0.02f
#define FLIRT_SCALE_MIN					+0.95f
#define FLIRT_SCALE_MAX					+1.05f
#define FLIRT_SKEW_STEP					+0.005f
#define FLIRT_SKEW_MIN					-0.05f
#define FLIRT_SKEW_MAX					+0.05f
#define FLIRT_SEARCH_PASS				2

#define FLIRT_HARDWARE_INTERPOLATION 1


//4mm stage
#define FLIRT_MINIMA_COUNT 3
#define FLIRT_HIGHER_RESOLUTION_START 1

#define FLIRT_BIN_WIDTH 1
#define FLIRT_BIN_COUNT (256/FLIRT_BIN_WIDTH)

#define FLIRT_BIN_WIDTH_4MM 1
#define FLIRT_BIN_COUNT_4MM (256/FLIRT_BIN_WIDTH_4MM)

#define FLIRT_BIN_WIDTH_8MM 2
#define FLIRT_BIN_COUNT_8MM (256/FLIRT_BIN_WIDTH_8MM)




#define WARP_SIZE 32
#define WARP_COUNT 8

#endif