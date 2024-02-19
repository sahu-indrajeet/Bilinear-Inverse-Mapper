//                                                                                                           //
//                                                                                                           //
//                                                                                                           //
//  ____  _                        _____                 _   _               _____             _             //
// / ___|| |__   __ _ _ __   ___  |  ___|   _ _ __   ___| |_(_) ___  _ __   | ____|_ __   __ _(_)_ __   ___  //
// \___ \| '_ \ / _` | '_ \ / _ \ | |_ | | | | '_ \ / __| __| |/ _ \| '_ \  |  _| | '_ \ / _` | | '_ \ / _ \ //
//  ___) | | | | (_| | |_) |  __/ |  _|| |_| | | | | (__| |_| | (_) | | | | | |___| | | | (_| | | | | |  __/ //
// |____/|_| |_|\__,_| .__/ \___| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_| |_____|_| |_|\__, |_|_| |_|\___| //
//                   |_|                                                                 |___/               //
//                                                                                                           //
//                                                                                                           //
//                                                                                                           //

/*
 * Library Name: Shape Function Engine
 * Author: Indrajeet Sahu
 * Version: 1.0
 * Date: February, 2024
 * Description:
                A powerful C library for finding the shape functions of any point with respect to a triangle or quadrilateral.
                Given a point and a polygon, the library runs to find the bilinear mapping parametric co-ordinates in case of quadrilateral and the Barycentric co-ordinates in case of a triangle.
 * License: MIT License
 *
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "shape_function_engine.h"

// todo: check if the code is up to date with the latest changes in your code in DEST

bool F_check_if_polygon_convex_ccw(double *polygon, size_t sides_num)
{    
    bool isConvex = true;
    bool isCCW = true;
    double v1x, v1y, v2x, v2y;  //  for edges vectors

    /*.........Check for CW vs CCW..............................*/
    if (IS_INBIL_POLYGN_CCW_CHECK)
    {
        // Check by choosing the first three vertex of the polygon
        
        // edge vector from 2nd to 3rd vertex
        v1x = polygon[4] - polygon[2];
        v1y = polygon[5] - polygon[3];

        // edge vector from 2nd to 1st vertex
        v2x = polygon[0] - polygon[2];
        v2y = polygon[1] - polygon[3];

        if ((v1x * v2y - v2x * v1y) <= 0)  // CCW check
        {
            isCCW = false;
        }
    }

    size_t array_size = 2 * sides_num;
    for (size_t i = 0; i < sides_num; i++)
    {
        // edge vector to the next point from this point
        v1x = polygon[(2 * i + 2) % array_size] - polygon[2 * i]; 
        v1y = polygon[(2 * i + 3) % array_size] - polygon[2 * i + 1];

        // edge vector to the previous point from this point
        v2x = polygon[(2 * i + array_size-2) % array_size] - polygon[2 * i];
        v2y = polygon[(2 * i + array_size-1) % array_size] - polygon[2 * i + 1];


        if( ( (v1x * v2x + v1y * v2y) / ( sqrt(v1x * v1x + v1y * v1y) * sqrt(v2x * v2x + v2y * v2y) ) ) < IS_INBIL_COS_MAX_ANGLE) //  cosine of angle between these two edges - Convexity check in a stricter sense than theoretical 180 degree
        {
            isConvex = false; 
            return false;
        }
    }
    return isConvex * isCCW ; 

}

void F_interpolate_quad(double QuadrilateralCoo[8], double u, double v, double *interpolated_x, double *interpolated_y)
{
    // Input Parametric Coordinates should be varying from 0 to 1 in both u and v
    *interpolated_x = (QuadrilateralCoo[0] * (1 - u) * (1 - v) +
                       QuadrilateralCoo[2] * (u)     * (1 - v) +
                       QuadrilateralCoo[4] * (u)     * (v)     +
                       QuadrilateralCoo[6] * (1 - u) * (v))    ;
    *interpolated_y = (QuadrilateralCoo[1] * (1 - u) * (1 - v) +
                       QuadrilateralCoo[3] * (u)     * (1 - v) +
                       QuadrilateralCoo[5] * (u)     * (v)     +
                       QuadrilateralCoo[7] * (1 - u) * (v))    ;
}


void F_InvBilMap_is_this_xi_trial_closer(double QuadrilateralCoo[8], double *xi1_new, double *xi2_new, double xi1_trial, double xi2_trial, double x, double y, double *current_difference)
{
    // Check if the (xi1_trial, xi2_trial) pair of parametric co-ordinates is closer to the exact solution. If so, use this as pair as the best solution so far.

    // if both xi1_trial and xi2_trial are valid quantities
    if (!isnan(xi1_trial) && !isnan(xi2_trial) &&
        !isinf(xi1_trial) && !isinf(xi2_trial)   )
    {        
        double interpolated_x, interpolated_y;
        F_interpolate_quad(QuadrilateralCoo, xi1_trial, xi2_trial, &interpolated_x, &interpolated_y);
        double trial_difference = (x - interpolated_x) * (x - interpolated_x) + (y - interpolated_y) * (y - interpolated_y);
        // if these parameters are giving lesser error than the previous ones then use these
        if (trial_difference < *current_difference)
        {
            *xi1_new = xi1_trial;
            *xi2_new = xi2_trial;
            *current_difference = trial_difference;
        }
    }
}

/*............ Inverse Mapping Function  ..............*/
bool F_InvBilMapQuadrilateral(double QuadrilateralCoo[8], double PointCoo[2], double* param_u, double* param_v)
{
    /* First check whether the quadrilateral is CCW and Convex */ 
    if (!F_check_if_polygon_convex_ccw(QuadrilateralCoo, 4))
    {
        // fprintf(stderr, "Unsuitable Quad - it's either CW or non-convex. Returning back..\n"); 
        return false;
    }

    /*  Parametric Space
    * 
    *   v                           v
    *   ^                where,     ^
    *   |                           |       
    *   4-------3                 +1|--------
    *   |       |                   |       |
    *   |       |                   |       |
    *   1-------2--->u              |---------->u
    *  O                           O       +1
    * 
    */

    /* Once confirmed that it is a suitable quad, then proceed further to calculate the bilinear mapping parameters*/
    
    double ISTOLER = IS_INBIL_ZERO_TOLER; // tolerance for near zero value comparison
    double ISQUADRTOL = IS_INBIL_QUADR_TOLER; // tolerance for Quadratic Equations which strictly check for point lying between [0,1] whereas some points might be slightly outside

    double x = PointCoo[0], y = PointCoo[1];
    size_t case_of_solution = 0;
    /*  Indices of CCW Quadrilateral
            x1, y1 : 0, 1
            x2, y2 : 2, 3
            x3, y3 : 4, 5
            x4, y4 : 6, 7           
    */

    /*... Axis aligned bounding box of the quadrilateral ...*/
    double max_x = fmax(fmax(QuadrilateralCoo[0], QuadrilateralCoo[2]),
                         fmax(QuadrilateralCoo[4], QuadrilateralCoo[6]));

    double min_x = fmin(fmin(QuadrilateralCoo[0], QuadrilateralCoo[2]),
                          fmin(QuadrilateralCoo[4], QuadrilateralCoo[6]));

    double max_y = fmax(fmax(QuadrilateralCoo[1], QuadrilateralCoo[3]),
                          fmax(QuadrilateralCoo[5], QuadrilateralCoo[7]));

    double min_y = fmin(fmin(QuadrilateralCoo[1], QuadrilateralCoo[3]),
                          fmin(QuadrilateralCoo[5], QuadrilateralCoo[7]));

    double diff_x = max_x - min_x;
    double diff_y = max_y - min_y;
    double min_dist = fmin(diff_x, diff_y);  // minimum side of the bounding box
    double allowable_error_distance = IS_INVBIL_DISTANCE_TOLER * min_dist; // physical distance of the allowable error

    double xi1, xi2;
    double    a1 = QuadrilateralCoo[0]                                                                    ;
    double    a2 = -QuadrilateralCoo[0] + QuadrilateralCoo[2]                                             ;
    double    a3 = -QuadrilateralCoo[0] + QuadrilateralCoo[6]                                             ;
    double    a4 = QuadrilateralCoo[0] - QuadrilateralCoo[2] + QuadrilateralCoo[4] - QuadrilateralCoo[6]  ;
    double    b1 = QuadrilateralCoo[1]                                                                    ;
    double    b2 = -QuadrilateralCoo[1] + QuadrilateralCoo[3]                                             ;
    double    b3 = -QuadrilateralCoo[1] + QuadrilateralCoo[7]                                             ;
    double    b4 = QuadrilateralCoo[1] - QuadrilateralCoo[3] + QuadrilateralCoo[5] - QuadrilateralCoo[7]  ;

    size_t Quadratic_Equation_to_be_used = 0;


    if (fabs(a4) < ISTOLER)
    {
        if (fabs(b4) < ISTOLER)
        {
            // Case 1.0 : a4 = 0 , b4 = 0
            xi1 = (x * b3 - y * a3 - a1 * b3 + b1 * a3) / (a2 * b3 - a3 * b2);
            xi2 = (x * b2 - y * a2 - a1 * b2 + a2 * b1) / (a3 * b2 - a2 * b3);          
            case_of_solution = 10;
        }
        else
        {
            // Case 2.0 : a4 = 0, b4 != 0
            if (fabs(a2) < ISTOLER)
            {
                if (fabs(a3) > ISTOLER)
                {
                    // Case 2.1 : a2 = 0 , a3 != 0
                    xi2 = (x - a1) / a3;
                    xi1 = (y - b1 - b3 * xi2) / (b2 + b4 * xi2);
                    case_of_solution = 21;
                }
            }
            else
            {
                if (fabs(a3) < ISTOLER)
                {
                    // Case 2.2 : a2 != 0, a3 = 0
                    xi1 = (x - a1) / a2;
                    xi2 = (y - b1 - b2 * xi1) / (b3 + b4 * xi1);
                    case_of_solution = 22;
                }
                else
                {
                    // Case 2.3 : a2 != 0, a3 != 0
                    if (fabs(a2) > fabs(a3))
                    {
                        Quadratic_Equation_to_be_used = 2;
                        case_of_solution = 23;
                    }
                    else
                    {
                        Quadratic_Equation_to_be_used = 1;
                        case_of_solution = 23;
                    }
                }
            }
            
        }
    }
    else
    {
        if (fabs(b4) < ISTOLER)
        {
            // Case 3.0 : a4 != 0, b4 = 0
            if (fabs(b2) < ISTOLER)
            {
                if (fabs(b3) > ISTOLER)
                {
                    // Case 3.1 : b2 = 0, b3 != 0
                    xi2 = (y - b1) / (b3);
                    xi1 = (x - a1 - a3 * xi2) / (a2 + a4 * xi2);
                    case_of_solution = 31;
                }
            }
            else
            {
                if (fabs(b3) < ISTOLER)
                {
                    // Case 3.2 : b2 != 0, b3 = 0
                    xi1 = (y - b1) / (b2);
                    xi2 = (x - a1 - a2 * xi1) / (a3 + a4 * xi1);
                    case_of_solution = 32;
                }
                else
                {
                    // Case 3.3 : b2 != 0, b3 != 0
                    if (fabs(b2) > fabs(b3))
                    {
                        Quadratic_Equation_to_be_used = 2;
                    }
                    else
                    {
                        Quadratic_Equation_to_be_used = 1;
                    }
                    case_of_solution = 33;
                }
            }
        }
        else
        {
            // Case 4.0 : a4 != 0, b4 != 0
            if (fabs(a2 * b4 - a4 * b2) < ISTOLER || fabs(a3 * b4 - a4 * b3) < ISTOLER)
            {
                if (fabs(a2 * b4 - a4 * b2) < fabs(a3 * b4 - a4 * b3))
                {
                    // Case 4.1 : a2 * b4 - a4 * b2 = 0
                    xi2 = (x * b4 - y * a4 - a1 * b4 + b1 * a4) / (a3 * b4 - a4 * b3);
                    if (fabs(a2 + a4 * xi2) > fabs(b2 + b4 * xi2))
                    {
                        xi1 = (x - a1 - a3 * xi2) / (a2 + a4 * xi2);
                    }
                    else
                    {
                        xi1 = (y - b1 - b3 * xi2) / (b2 + b4 * xi2);
                    }
                    case_of_solution = 41;
                }
                else
                {
                    // Case 4.2 : a3 * b4 - a4 * b3 = 0
                    xi1 = (x * b4 - y * a4 - a1 * b4 + a4 * b1) / (a2 * b4 - a4 * b2);
                    if (fabs(a3 + a4 * xi1) > fabs(b3 + b4 * xi1))
                    {
                        xi2 = (x - a1 - a2 * xi1) / (a3 + a4 * xi1);
                    }
                    else
                    {
                        xi2 = (y - b1 - b2 * xi1) / (b3 + b4 * xi1);
                    }
                    case_of_solution = 42;
                }
            }
            else
            {
                // Case 4.3 :  a2 * b4 - a4 * b2 != 0 and  a3 * b4 - a4 * b3 != 0
                if (fabs(a3 * b4 - a4 * b3) > fabs(a2 * b4 - a4 * b2))
                {
                    Quadratic_Equation_to_be_used = 1;
                }
                else
                {
                    Quadratic_Equation_to_be_used = 2;
                }
                case_of_solution = 43;
            }            
        }
    }

    if (Quadratic_Equation_to_be_used == 1)
    {
        double P = a4 * b2 - a2 * b4;
        double Q = b4 * x - a4 * y + a4 * b1 - a1 * b4 + a3 * b2 - a2 * b3;
        double R = b3 * x - a3 * y + a3 * b1 - a1 * b3;

        double discrim_root = sqrt(Q * Q - 4 * P * R);

        xi1 = (-Q + discrim_root) / (2 * P);

        if (xi1 < -ISQUADRTOL || xi1 >(1 + ISQUADRTOL))
        {
            xi1 = (-Q - discrim_root) / (2 * P);
        }

        if (fabs(a3 + a4 * xi1) > fabs(b3 + b4 * xi1))
        {
            xi2 = (x - a1 - a2 * xi1) / (a3 + a4 * xi1);
        }
        else
        {
            xi2 = (y - b1 - b2 * xi1) / (b3 + b4 * xi1);
        }
    }

    if (Quadratic_Equation_to_be_used == 2)
    {
        double S = a4 * b3 - a3 * b4;
        double T = b4 * x - a4 * y + a4 * b1 - a1 * b4 - a3 * b2 + a2 * b3;
        double U = b2 * x - a2 * y + a2 * b1 - a1 * b2;

        double discrim_root = sqrt(T * T - 4 * S * U);

        xi2 = (-T + discrim_root) / (2 * S);

        if (xi2 < -ISQUADRTOL || xi2 >(1 + ISQUADRTOL))
        {
            xi2 = (-T - discrim_root) / (2 * S);
        }

        if (fabs(a2 + a4 * xi2) > fabs(b2 + b4 * xi2))
        {
            xi1 = (x - a1 - a3 * xi2) / (a2 + a4 * xi2);
        }
        else
        {
            xi1 = (y - b1 - b3 * xi2) / (b2 + b4 * xi2);
        }
    }


    double interpolated_x = (QuadrilateralCoo[0] * (1 - xi1) * (1 - xi2) +
                                QuadrilateralCoo[2] * (xi1) * (1 - xi2) +
                                QuadrilateralCoo[4] * (xi1) * (xi2)+
                                QuadrilateralCoo[6] * (1 - xi1) * (xi2));
    double interpolated_y = (QuadrilateralCoo[1] * (1 - xi1) * (1 - xi2) +
                                QuadrilateralCoo[3] * (xi1) * (1 - xi2) +
                                QuadrilateralCoo[5] * (xi1) * (xi2)+
                                QuadrilateralCoo[7] * (1 - xi1) * (xi2));

    double difference = (x - interpolated_x) * (x - interpolated_x) + (y - interpolated_y) * (y - interpolated_y);

    // if the error is beyond the satisfactory limits then go and check all other cases and choose the one with least distance from the actual point as the solution
    if (difference > allowable_error_distance)
    {
        double xi1_new = xi1;
        double xi2_new = xi2;
        double xi1_trial;
        double xi2_trial;

        double current_difference = difference;

        // Case 1.0 : a4 = 0, b4 = 0
        if (case_of_solution != 10)
        {
            xi1_trial = (x * b3 - y * a3 - a1 * b3 + b1 * a3) / (a2 * b3 - a3 * b2);
            xi2_trial = (x * b2 - y * a2 - a1 * b2 + a2 * b1) / (a3 * b2 - a2 * b3);
            // if these parameters have errors less than the case of previous parameters then these new parameters would be used
            F_InvBilMap_is_this_xi_trial_closer(QuadrilateralCoo, &xi1_new, &xi2_new, xi1_trial, xi2_trial, x, y, &current_difference);
        }

        // Case 2.1 : a2 = 0 , a3 != 0
        if (case_of_solution != 21)
        {
            xi2_trial = (x - a1) / a3;
            xi1_trial = (y - b1 - b3 * xi2) / (b2 + b4 * xi2_trial);
            F_InvBilMap_is_this_xi_trial_closer(QuadrilateralCoo, &xi1_new, &xi2_new, xi1_trial, xi2_trial, x, y, &current_difference);
        }

        // Case 2.2 : a2 != 0, a3 = 0
        if (case_of_solution != 22)
        {
            xi1_trial = (x - a1) / a2;
            xi2_trial = (y - b1 - b2 * xi1_trial) / (b3 + b4 * xi1_trial);
            F_InvBilMap_is_this_xi_trial_closer(QuadrilateralCoo, &xi1_new, &xi2_new, xi1_trial, xi2_trial, x, y, &current_difference);
        }

        // Case 3.1 : b2 = 0, b3 != 0
        if (case_of_solution != 31)
        {
            xi2 = (y - b1) / (b3);
            xi1_trial = (x - a1 - a3 * xi2) / (a2 + a4 * xi2);
            F_InvBilMap_is_this_xi_trial_closer(QuadrilateralCoo, &xi1_new, &xi2_new, xi1_trial, xi2_trial, x, y, &current_difference);
        }

        // Case 3.2 : b2 != 0, b3 = 0
        if (case_of_solution != 32)
        {
            xi1_trial = (y - b1) / (b2);
            xi2_trial = (x - a1 - a2 * xi1_trial) / (a3 + a4 * xi1_trial);
            F_InvBilMap_is_this_xi_trial_closer(QuadrilateralCoo, &xi1_new, &xi2_new, xi1_trial, xi2_trial, x, y, &current_difference);
        }

        // Case 4.1 : a2 * b4 - a4 * b2 = 0
        if (case_of_solution != 41)
        {
            xi2_trial = (x * b4 - y * a4 - a1 * b4 + b1 * a4) / (a3 * b4 - a4 * b3);
            if (fabs(a2 + a4 * xi2_trial) > fabs(b2 + b4 * xi2_trial))
            {
                xi1_trial = (x - a1 - a3 * xi2_trial) / (a2 + a4 * xi2_trial);
            }
            else
            {
                xi1_trial = (y - b1 - b3 * xi2_trial) / (b2 + b4 * xi2_trial);
            }
            F_InvBilMap_is_this_xi_trial_closer(QuadrilateralCoo, &xi1_new, &xi2_new, xi1_trial, xi2_trial, x, y, &current_difference); 
        }

        // Case 4.2 : a3 * b4 - a4 * b3 = 0
        if (case_of_solution != 42)
        {
            xi1_trial = (x * b4 - y * a4 - a1 * b4 + a4 * b1) / (a2 * b4 - a4 * b2);
            if (fabs(a3 + a4 * xi1_trial) > fabs(b3 + b4 * xi1_trial))
            {
                xi2_trial = (x - a1 - a2 * xi1_trial) / (a3 + a4 * xi1_trial);
            }
            else
            {
                xi2_trial = (y - b1 - b2 * xi1_trial) / (b3 + b4 * xi1_trial);
            }
            F_InvBilMap_is_this_xi_trial_closer(QuadrilateralCoo, &xi1_new, &xi2_new, xi1_trial, xi2_trial, x, y, &current_difference); 
        }

        // First Quadratic 
        if (Quadratic_Equation_to_be_used != 1)
        {
            double P = a4 * b2 - a2 * b4;
            double Q = b4 * x - a4 * y + a4 * b1 - a1 * b4 + a3 * b2 - a2 * b3;
            double R = b3 * x - a3 * y + a3 * b1 - a1 * b3;

            double discrim_root = sqrt(Q * Q - 4 * P * R);
            xi1_trial = (-Q + discrim_root) / (2 * P);

            if (xi1_trial < -ISQUADRTOL || xi1_trial >(1 + ISQUADRTOL))
            {
                xi1_trial = (-Q - discrim_root) / (2 * P);
            }

            if (fabs(a3 + a4 * xi1_trial) > fabs(b3 + b4 * xi1_trial))
            {
                xi2_trial = (x - a1 - a2 * xi1_trial) / (a3 + a4 * xi1_trial);
            }
            else
            {
                xi2_trial = (y - b1 - b2 * xi1_trial) / (b3 + b4 * xi1_trial);
            }
            F_InvBilMap_is_this_xi_trial_closer(QuadrilateralCoo, &xi1_new, &xi2_new, xi1_trial, xi2_trial, x, y, &current_difference);
        }

        // Second Quadratic
        if (Quadratic_Equation_to_be_used != 2)
        {
            double S = a4 * b3 - a3 * b4;
            double T = b4 * x - a4 * y + a4 * b1 - a1 * b4 - a3 * b2 + a2 * b3;
            double U = b2 * x - a2 * y + a2 * b1 - a1 * b2;

            double discrim_root = sqrt(T * T - 4 * S * U);

            xi2_trial = (-T + discrim_root) / (2 * S);

            if (xi2_trial < -ISQUADRTOL || xi2_trial >(1 + ISQUADRTOL))
            {
                xi2_trial = (-T - discrim_root) / (2 * S);
            }

            if (fabs(a2 + a4 * xi2_trial) > fabs(b2 + b4 * xi2_trial))
            {
                xi1_trial = (x - a1 - a3 * xi2_trial) / (a2 + a4 * xi2_trial);
            }
            else
            {
                xi1_trial = (y - b1 - b3 * xi2_trial) / (b2 + b4 * xi2_trial);
            }
            F_InvBilMap_is_this_xi_trial_closer(QuadrilateralCoo, &xi1_new, &xi2_new, xi1_trial, xi2_trial, x, y, &current_difference);
        }

        // After having gone through all the cases, use the best parameters that give least error
        xi1 = xi1_new;
        xi2 = xi2_new;

    }

    if ( isnan(xi1) || isinf(xi1) ||
         isnan(xi2) || isinf(xi2)        )
    {
        fprintf(stderr, "F_InvBilMapQuadrilateral : isnan or inf error. Returning back...\n"); // displays on the running console 
        return false; // shouldn't ever come here if it is non-convex quad
    }

    if (IS_INBIL_POINT_LEEWAY > 0.0)
    {
        // if xi1 or xi2 are just outside the parametric domain then make them to be on the edge of it. It can possibly happen due to some numerical precision issues or the point that came to this function was actually outside the quad
        double allowable_min = IS_INBIL_POINT_LEEWAY;
        double allowable_max = 1 + IS_INBIL_POINT_LEEWAY;
        if (xi1 < 0.0 && xi1 > allowable_min)
        {
            xi1 = 0.0;
        }
        if (xi1 > 1.0 && xi1 < allowable_max)
        {
            xi1 = 1.0;
        }
        if (xi2 < 0.0 && xi2 > allowable_min)
        {
            xi2 = 0.0;
        }
        if (xi2 > 1.0 && xi2 < allowable_max)
        {
            xi2 = 1.0;
        }
    }

    // Finaly, we have got our parameter values that this function was intended for
    *param_u = xi1;
    *param_v = xi2;

    return true;
}

// to convert the non-standard parametric coordinates [0,1] to the standard parametric coordinates [-1,1]
void F_Quad_Std_ParamCoo(double u, double v, double* xi1, double* xi2)
{
    /*
    * both u   and v    belong to [ 0,1]
    * both xi1 and xi2  belong to [-1,1]
    * 
    *   Non-standard Parametric Space
    *
    *   v                           v
    *   ^                where,     ^
    *   |                           |       
    *   4-------3                 +1|--------
    *   |       |                   |       |
    *   |       |                   |       |
    *   1-------2--->u              |---------->u
    *  O                           O       +1
    * 
    * 
    *  Standard Parametric Space
    * 
    *         xi2
    *         /\
    *         |
    *         |
    * 4---------------3
    * |       |       |
    * |       |       |
    * |-------o-------|---->xi1     where, xi1 and xi2 \in  [-1,1]
    * |       |       |
    * |       |       |
    * 1---------------2
    * 
    */

    *xi1 = 2 * u - 1;
    *xi2 = 2 * v - 1;
}

void F_ShapeFunction_values(double xi1, double xi2, double Psi[4])
{
   /*
    *         xi2
    *         /\
    *         |
    *         |
    * 4---------------3
    * |       |       |
    * |       |       |
    * |-------o-------|---->xi1
    * |       |       |
    * |       |       |
    * 1---------------2
    * 
    */
    // bilinear shape functions of the quadrilateral
    Psi[0] = 0.25 * (1 - xi1) * (1 - xi2);
    Psi[1] = 0.25 * (1 + xi1) * (1 - xi2);
    Psi[2] = 0.25 * (1 + xi1) * (1 + xi2);
    Psi[3] = 0.25 * (1 - xi1) * (1 + xi2);

}

void Shape_function_engine_Quad(double Quadrilateral[8], double PointCoo[2], double Shape_func_values[4])
{
	// Incoming Quadrilateral must be CCW 
	double temp_u, temp_v;
	double xi, eta;
    F_InvBilMapQuadrilateral(Quadrilateral, PointCoo, &temp_u, &temp_v);
    F_Quad_Std_ParamCoo(temp_u, temp_v, &xi, &eta);	
    F_ShapeFunction_values(xi, eta, Shape_func_values);	

}

// area of a triangle (assuming CCW ordering)
double F_triangle_area(double Triangle[6])
{
    double x1, x2, x3;
    double y1, y2, y3;
    x1 = Triangle[0];
    x2 = Triangle[2];
    x3 = Triangle[4];
    y1 = Triangle[1];
    y2 = Triangle[3];
    y3 = Triangle[5];

    return 0.5 * ((x1 * y2 + x2 * y3 + x3 * y1) - (x2 * y1 + x3 * y2 + x1 * y3));  // will be positive for CCW triangle, and negative for CW triangle
} 

// To find the Barycentric Co-oridnates of the given point wrt the triangle
bool F_Triangle_Pt_BCC(double Triangle[6], double PointCoo[2], double *u, double *v, double *w)
{
    /* First check whether the triangle is CCW and Convex */
    if (!F_check_if_polygon_convex_ccw(Triangle, 3))
    {
        // fprintf(stderr, "Unsuitable Triangle - it's either CW or non-convex. Returning back..\n"); 
        return false;
    }

    double Area_subtriangles[3]; // area of the three sub-triangles formed by the given triangle and the given point
    double Area_wholetriangle = 0;          //area of the whole triangle

    double subtriangle[6];
    subtriangle[0] = PointCoo[0];
    subtriangle[1] = PointCoo[1];

    // loop over all sub-triangles
    for (size_t iter = 1; iter < 4; iter++)
    {
        subtriangle[2] = Triangle[(2 * iter) % 6];
        subtriangle[3] = Triangle[(2 * iter + 1) % 6];

        subtriangle[4] = Triangle[(2 * iter + 2) % 6];
        subtriangle[5] = Triangle[(2 * iter + 1 + 2) % 6];

        // area of this sub-triangle
        Area_subtriangles[iter - 1] = F_triangle_area(subtriangle);
        Area_wholetriangle += Area_subtriangles[iter - 1];
    }

    *u = Area_subtriangles[0] / Area_wholetriangle;
    *v = Area_subtriangles[1] / Area_wholetriangle;
    *w = Area_subtriangles[2] / Area_wholetriangle;
}
