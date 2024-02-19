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

/* ............ User Controls ............*/

// polygon rejection if it's not CCW ordered, or it's not sufficiently convex
#define IS_INBIL_POLYGN_CCW_CHECK   true
#define IS_INBIL_COS_MAX_ANGLE      -0.98481    // cosine of maximum angle between any two sides of the polygon for convex check. Program returns false if any polygon corner exceeds this angle

// quadrilateral inverse mapping
#define IS_INVBIL_DISTANCE_TOLER    1e-20       // how much maximum part of minimum distance in quad do you want to allow as error

#define IS_INBIL_ZERO_TOLER         1e-12       // tolerance for zero value comparison
#define IS_INBIL_QUADR_TOLER        1e-5        // tolerance for use in roots of quadratic equation

#define IS_INBIL_POINT_LEEWAY       0.1         // refers to the parametric distance in both directions (\Delta u or \Delta v) for the permitted point's distance outside the quadrilateral. If the point lies within this parametric distance, the inverse mapping function would return its parametric coordinates to be lying on the edge.

/*....................................*/



    /*  Indices of CCW Quadrilateral
            x1, y1 : 0, 1
            x2, y2 : 2, 3
            x3, y3 : 4, 5
            x4, y4 : 6, 7
        Indices for PointCoo
            x, y : 0, 1
    */

// Function to be used to find the point of interpolation of quadrilateral using non-standard parametric co-ordinates (u,v)
extern void F_interpolate_quad(double QuadrilateralCoo[8], double u, double v, double* interpolated_x, double* interpolated_y);

// Function to be used to calculate the parametric co-ordinates of bilinear mapping with respect to a quadrilateral
extern bool F_InvBilMapQuadrilateral(double QuadrilateralCoo[8], double PointCoo[2], double* param_u, double* param_v);

// Function to be used to calculate the shape functions
extern void Shape_function_engine_Quad(double Quadrilateral[8], double PointCoo[2], double Shape_func_values[4]);


    /*  Indices of CCW Triangle
            x1, y1 : 0, 1
            x2, y2 : 2, 3
            x3, y3 : 4, 5
        Indices for PointCoo
            x, y : 0, 1
    */

extern bool F_Triangle_Pt_BCC(double Triangle[6], double PointCoo[2], double* u, double* v, double* w);

