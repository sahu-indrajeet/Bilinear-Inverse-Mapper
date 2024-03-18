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

#include <stdbool.h>
#include <math.h>

#include "bilinear_inverse_mapper.h"


bool F_check_if_polygon_convex(double *polygon, size_t sides_num)
{   
    // Assuming a CCW ordered quadrilateral, it returns whether the given quadrilateral is convex(true) or not(false)

    // convexity check requires that 1) cross product of   next x previous   edge vectors at all corners are +ve
    // and, for more strict convexity check, the (cosine of) angle at all corners should be less than the (cosine of) maximum angle

    double v1x, v1y, v2x, v2y;  //  for edges vectors

    size_t array_size = 2 * sides_num;
    for (size_t i = 0; i < sides_num; i++)
    {
        // edge vector to the next point from this point
        v1x = polygon[(2 * i + 2) % array_size] - polygon[2 * i]; 
        v1y = polygon[(2 * i + 3) % array_size] - polygon[2 * i + 1];

        // edge vector to the previous point from this point
        v2x = polygon[(2 * i + array_size-2) % array_size] - polygon[2 * i];
        v2y = polygon[(2 * i + array_size-1) % array_size] - polygon[2 * i + 1];

        // next x previous cross-product
        if ((v1x * v2y - v2x * v1y) <= 0)  //if any cross-product is going into the screen
        {
            return false; // this vertex is bringing non-convexity to the quadrilateral
        }
        if( ( (v1x * v2x + v1y * v2y) / ( sqrt(v1x * v1x + v1y * v1y) * sqrt(v2x * v2x + v2y * v2y) ) ) < IS_INBIL_COS_MAX_ANGLE) //  cosine of angle between these two edges - Convexity check in a stricter sense than theoretical 180 degree
        {
            return false;
        }
    }
    // it's convex with all angles being less than the threshold
    return true;
}


bool F_check_if_polygon_min_angle(double* polygon, size_t sides_num)
{
    // just checks for the minimum internal or external angle cosine
    // returns true if ok but returns false if not ok (e.g. degenerate quad)

    double v1x, v1y, v2x, v2y;  //  for edges vectors

    size_t array_size = 2 * sides_num;
    for (size_t i = 0; i < sides_num; i++)
    {
        // edge vector to the next point from this point
        v1x = polygon[(2 * i + 2) % array_size] - polygon[2 * i];
        v1y = polygon[(2 * i + 3) % array_size] - polygon[2 * i + 1];

        // edge vector to the previous point from this point
        v2x = polygon[(2 * i + array_size - 2) % array_size] - polygon[2 * i];
        v2y = polygon[(2 * i + array_size - 1) % array_size] - polygon[2 * i + 1];

        if (((v1x * v2x + v1y * v2y) / (sqrt(v1x * v1x + v1y * v1y) * sqrt(v2x * v2x + v2y * v2y))) < IS_INBIL_COS_MAX_ANGLE) //  cosine of angle between these two edges - Convexity check in a stricter sense than theoretical 180 degree
        {
            return false;
        }
    }
    // it's convex with all angles being less than the threshold
    return true;
}

void F_interpolate_quad(double QuadrilateralCoo[8], double xi1, double xi2, double *interpolated_x, double *interpolated_y)
{
    // Input Parametric Coordinates should be varying from 0 to 1 in both u and v
    *interpolated_x = (QuadrilateralCoo[0] * (1 - xi1) * (1 - xi2) +
                       QuadrilateralCoo[2] * (xi1)     * (1 - xi2) +
                       QuadrilateralCoo[4] * (xi1)     * (xi2)     +
                       QuadrilateralCoo[6] * (1 - xi1) * (xi2))    ;
    *interpolated_y = (QuadrilateralCoo[1] * (1 - xi1) * (1 - xi2) +
                       QuadrilateralCoo[3] * (xi1)     * (1 - xi2) +
                       QuadrilateralCoo[5] * (xi1)     * (xi2)     +
                       QuadrilateralCoo[7] * (1 - xi1) * (xi2))    ;
}



void F_quadr_roots(long double a,long double b, double c,long double roots[2])
{
    // function to calculate the roots of a quadratic equation

    long double discrim = b * b - 4 * a * c;
    if (b > 0)
    {
        roots[0] = (-2 * c) / (b + sqrt(discrim)); // smaller root
        roots[1] = (-b - sqrt(discrim)) / (2*a);
        return;
    }
    roots[0] = (2 * c) / (-b + sqrt(discrim)); // smaller root
    roots[1] = (b - sqrt(discrim)) / (-2 * a);
}


/*............ Inverse Mapping Function  ..............*/
bool F_InBilMap_Quadrilateral(double QuadrilateralCoo[8], double PointCoo[2], double* param_xi1, double* param_xi2, bool *both_roots_wanted)
{

    /* First check whether the quadrilateral is CCW and Convex */ 
#if IS_INBIL_POLYGN_CONVEX_CHECK
    if (!F_check_if_polygon_convex(QuadrilateralCoo, 4))
    {
        // fprintf(stderr, "Unsuitable Quad - it's either CW or non-convex. Returning back..\n"); 
        return false;
    }
    /* Once confirmed that it is a suitable quad, then proceed further to calculate the bilinear mapping parameters*/
#endif

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

    double x = PointCoo[0], y = PointCoo[1];
    /*  Indices of CCW Quadrilateral
            x1, y1 : 0, 1
            x2, y2 : 2, 3
            x3, y3 : 4, 5
            x4, y4 : 6, 7           
    */


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
    double P, S;
    P = a4 * b2 - a2 * b4;
    S = a4 * b3 - a3 * b4 ;

    // Case 4.0 : a4 != 0, b4 != 0
    if (fabs(P) > IS_INBIL_ZERO_TOLER || fabs(S) > IS_INBIL_ZERO_TOLER)  
    {
        // Case 4.3 :  a2 * b4 - a4 * b2 != 0 and  a3 * b4 - a4 * b3 != 0
        if (fabs(S) < fabs(P))
        {
            Quadratic_Equation_to_be_used = 1;
        }
        else
        {
            Quadratic_Equation_to_be_used = 2;
        }        
    }
    else
    {
        // both a4 and b4 are non-zero
        if (fabs(a4) > IS_INBIL_ZERO_TOLER && fabs(b4) > IS_INBIL_ZERO_TOLER)
        {
            if (fabs(P) < fabs(S))
            {
                // Case 4.1 : a2 * b4 - a4 * b2 = 0
                xi2 = (x * b4 - y * a4 - a1 * b4 + b1 * a4) / (-S);
                if (fabs(a2 + a4 * xi2) > fabs(b2 + b4 * xi2))
                {
                    xi1 = (x - a1 - a3 * xi2) / (a2 + a4 * xi2);
                }
                else
                {
                    xi1 = (y - b1 - b3 * xi2) / (b2 + b4 * xi2);
                }
            }
            else
            {
                // Case 4.2 : a3 * b4 - a4 * b3 = 0
                xi1 = (x * b4 - y * a4 - a1 * b4 + a4 * b1) / (-P);
                if (fabs(a3 + a4 * xi1) > fabs(b3 + b4 * xi1))
                {
                    xi2 = (x - a1 - a2 * xi1) / (a3 + a4 * xi1);
                }
                else
                {
                    xi2 = (y - b1 - b2 * xi1) / (b3 + b4 * xi1);
                }
            }
        }
        else if (fabs(a4) > IS_INBIL_ZERO_TOLER && fabs(b4) < IS_INBIL_ZERO_TOLER)
        {
            // Case 3.0 : a4 != 0, b4 = 0
            if (fabs(b3) > fabs(b2))
            {
                // Case 3.1 : b2 = 0, b3 != 0
                xi2 = (y - b1) / (b3);
                xi1 = (x - a1 - a3 * xi2) / (a2 + a4 * xi2);
            }
            else
            {
                // Case 3.2 : b2 != 0, b3 = 0
                xi1 = (y - b1) / (b2);
                xi2 = (x - a1 - a2 * xi1) / (a3 + a4 * xi1);
            }
        }
        else if (fabs(a4) < IS_INBIL_ZERO_TOLER && fabs(b4) > IS_INBIL_ZERO_TOLER)
        {
            // Case 2.0 : a4 = 0, b4 != 0
            if (fabs(a2) < fabs(a3))
            {
                // Case 2.1 : a2 = 0 , a3 != 0
                xi2 = (x - a1) / a3;
                xi1 = (y - b1 - b3 * xi2) / (b2 + b4 * xi2);
            }
            else
            {
                // Case 2.2 : a2 != 0, a3 = 0
                xi1 = (x - a1) / a2;
                xi2 = (y - b1 - b2 * xi1) / (b3 + b4 * xi1);
            }
        }
        else 
        {
            // Case 1.0 : a4 = 0 , b4 = 0
            xi1 = (x * b3 - y * a3 - a1 * b3 + b1 * a3) / (a2 * b3 - a3 * b2);
            xi2 = (x * b2 - y * a2 - a1 * b2 + a2 * b1) / (a3 * b2 - a2 * b3);
        }
    }


    if (Quadratic_Equation_to_be_used == 1)
    {
        double Q = b4 * x - a4 * y + a4 * b1 - a1 * b4 + a3 * b2 - a2 * b3;
        double R = b3 * x - a3 * y + a3 * b1 - a1 * b3;

        double roots[2];
        F_quadr_roots(P, Q, R, roots);

        if (*both_roots_wanted)
        {
            *param_xi1       = roots[0];
            *(param_xi1 + 1) = roots[1];

            xi1 = roots[0];
            if (fabs(a3 + a4 * xi1) > fabs(b3 + b4 * xi1))
            {
                *param_xi2       = (x - a1 - a2 * xi1) / (a3 + a4 * xi1);
            }
            else
            {
                *param_xi2 = (y - b1 - b2 * xi1) / (b3 + b4 * xi1);
            }
            xi1 = roots[1];
            if (fabs(a3 + a4 * xi1) > fabs(b3 + b4 * xi1))
            {
                *(param_xi2 + 1) = (x - a1 - a2 * xi1) / (a3 + a4 * xi1);
            }
            else
            {
                *(param_xi2 + 1) = (y - b1 - b2 * xi1) / (b3 + b4 * xi1);
            }
            //if (isnan(*param_u) || isinf(*(param_u + 1)) ||
            //    isnan(*param_v) || isinf(*(param_v + 1)))
            //{
            //    fprintf(stderr, "F_InvBilMapQuadrilateral : isnan or inf error. Returning back...\n"); // displays on the running console 
            //    return false; // shouldn't ever come here if it is non-convex quad
            //}
           return true;
        }
        if (roots[0] < -IS_INBIL_QUADR_TOLER || roots[0] > (1 + IS_INBIL_QUADR_TOLER))
        {
            xi1 = roots[1];
        }
        else
        {
            xi1 = roots[0];
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
        double T = b4 * x - a4 * y + a4 * b1 - a1 * b4 - a3 * b2 + a2 * b3;
        double U = b2 * x - a2 * y + a2 * b1 - a1 * b2;

        double roots[2];
        F_quadr_roots(S, T, U, roots);

        if (*both_roots_wanted)
        {
            *param_xi2       = roots[0];
            *(param_xi2 + 1) = roots[1];

            xi2 = roots[0];
            if (fabs(a2 + a4 * xi2) > fabs(b2 + b4 * xi2))
            {
                *param_xi1       = (x - a1 - a3 * xi2) / (a2 + a4 * xi2);
            }
            else
            {
                *param_xi1 = (y - b1 - b3 * xi2) / (b2 + b4 * xi2);
            }
            xi2 = roots[1];
            if (fabs(a2 + a4 * xi2) > fabs(b2 + b4 * xi2))
            {
                *(param_xi1 + 1) = (x - a1 - a3 * xi2) / (a2 + a4 * xi2);
            }
            else
            {
                *(param_xi1 + 1) = (y - b1 - b3 * xi2) / (b2 + b4 * xi2);
            }
            //if (isnan(*param_u) || isinf(*(param_u+1)) ||
            //    isnan(*param_v) || isinf(*(param_v+1)) )
            //{
            //    fprintf(stderr, "F_InvBilMapQuadrilateral : isnan or inf error. Returning back...\n"); // displays on the running console 
            //    return false; // shouldn't ever come here if it is non-convex quad
            //}
            return true;
        }

        if (roots[0] < -IS_INBIL_QUADR_TOLER || roots[0] > (1 + IS_INBIL_QUADR_TOLER))
        {
            xi2 = roots[1];
        }
        else
        {
            xi2 = roots[0];
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

    //if (isnan(xi1) || isinf(xi1) ||
    //    isnan(xi2) || isinf(xi2))
    //{
    //    fprintf(stderr, "F_InvBilMapQuadrilateral : isnan or inf error. Returning back...\n"); // displays on the running console 
    //    return false; // shouldn't ever come here if it is non-convex quad
    //}

    *param_xi1 = xi1;
    *param_xi2 = xi2;
    *both_roots_wanted = false;
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

void F_InterpolFunction_values(double xi1, double xi2, double Psi[4])
{
   /*
    *         xi2
    *         /\
    *         |
    *         |
    *     +1  4-------3
    *         |       |
    *         |       |
    *         1-------2---->xi1
    *       O         +1
    * 
    */
    // bilinear shape functions of the quadrilateral
    Psi[0] =  (1 - xi1) * (1 - xi2);
    Psi[1] =  (    xi1) * (1 - xi2);
    Psi[2] =  (    xi1) * (    xi2);
    Psi[3] =  (1 - xi1) * (    xi2);

}

void F_Interpol_Function_Calculator_Quad(double Quadrilateral[8], double PointCoo[2], double Shape_func_values[4])
{
    // Calcualtes the values of interpolation functions at the given point internal to the quadrilateral (only suitable for single root if quadrilateral with non-parallel opposite sides)

	// Incoming Quadrilateral must be CCW 
	double xi1, xi2;
    bool both_roots_wanted = false;
    F_InBilMap_Quadrilateral(Quadrilateral, PointCoo, &xi1, &xi2, &both_roots_wanted);
    F_InterpolFunction_values(xi1, xi2, Shape_func_values);	
}

double F_triangle_area(double Triangle[6])
{
    // area of a triangle (assuming CCW ordering)

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

bool F_Triangle_Pt_BCC(double Triangle[6], double PointCoo[2], double *u, double *v, double *w)
{
    // To find the Barycentric Co-oridnates of the given point wrt the triangle
        
    /* First check whether the triangle is CCW and Convex */
    if (!F_check_if_polygon_convex(Triangle, 3))
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

    return true;
}

