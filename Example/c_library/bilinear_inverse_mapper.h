/**
* /======================================================================================================================\
* ||  ____  _ _ _                           ___                                    __  __                               ||
* || | __ )(_) (_)_ __   ___  __ _ _ __    |_ _|_ ____   _____ _ __ ___  ___      |  \/  | __ _ _ __  _ __   ___ _ __   ||
* || |  _ \| | | | '_ \ / _ \/ _` | '__|____| || '_ \ \ / / _ \ '__/ __|/ _ \_____| |\/| |/ _` | '_ \| '_ \ / _ \ '__|  ||
* || | |_) | | | | | | |  __/ (_| | | |_____| || | | \ V /  __/ |  \__ \  __/_____| |  | | (_| | |_) | |_) |  __/ |     ||
* || |____/|_|_|_|_| |_|\___|\__,_|_|      |___|_| |_|\_/ \___|_|  |___/\___|     |_|  |_|\__,_| .__/| .__/ \___|_|     ||
* ||                                                                                           |_|   |_|                ||
* \======================================================================================================================/
*
*
* Library Name: Bilinear-Inverse-Mapper
* Author: Indrajeet Sahu
* Version: 1.0
* Date: March, 2024
* Description:
               A powerful C library for finding the interpolation functions of any point with respect to a triangle or quadrilateral.
               Given a point and a polygon, the library runs to find the bilinear mapping parametric co-ordinates in case of quadrilateral and the Barycentric co-ordinates in case of a triangle.
* License: MIT License
*
*/

/* ............ User Controls ............*/

// polygon rejection if it's not CCW ordered, or it's not sufficiently convex
#define IS_INBIL_POLYGN_CONVEX_CHECK false       // to check if the polygon is convex or not (it assumes that the polygon is CCW ordered)
#define IS_INBIL_COS_MAX_ANGLE    -0.98481    // cosine of maximum angle between any two sides of the polygon for convex check. Program returns false if any polygon corner exceeds this angle

// tolerance for use in roots of quadratic equation to select desirable root between 0 and 1 only
#define IS_INBIL_QUADR_TOLER        1e-10        

/*....................................*/


#include <float.h>
#define IS_INBIL_ZERO_TOLER    DBL_EPSILON       // tolerance for zero value comparison



    /*  Indices of CCW Quadrilateral
            x1, y1 : 0, 1
            x2, y2 : 2, 3
            x3, y3 : 4, 5
            x4, y4 : 6, 7
        Indices for PointCoo
            x, y : 0, 1
    */

// Function to be used to calculate the parametric co-ordinates of bilinear mapping with respect to a quadrilateral
extern bool F_InBilMap_Quadrilateral(double QuadrilateralCoo[8], double PointCoo[2], double* param_xi1, double* param_xi2, bool* both_roots_wanted);

// Function to be used to calculate the shape functions
extern void F_Interpol_Function_Calculator_Quad(double Quadrilateral[8], double PointCoo[2], double Shape_func_values[4]);

// Function to be used to find the point of interpolation of quadrilateral using non-standard parametric co-ordinates (u,v) with each coordinate beloning to [0,1]
extern void F_interpolate_quad(double QuadrilateralCoo[8], double xi1, double xi2, double* interpolated_x, double* interpolated_y);



    /*  Indices of CCW Triangle
            x1, y1 : 0, 1
            x2, y2 : 2, 3
            x3, y3 : 4, 5
        Indices for PointCoo
            x, y : 0, 1
    */

extern bool F_Triangle_Pt_BCC(double Triangle[6], double PointCoo[2], double* u, double* v, double* w);

