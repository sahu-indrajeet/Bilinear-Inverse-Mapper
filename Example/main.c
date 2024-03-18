

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "bilinear_inverse_mapper.h"

int main()
{
    double QuadrilateralCoo[8];
    double PointCoo[2];
    double xi1, xi2; // parametric coordinates of the given point
    bool both_roots_wanted;

    /********************************************/
    /******* For only one root ***********/
    /********************************************/
    both_roots_wanted = false;

    QuadrilateralCoo[0] = -1;
    QuadrilateralCoo[1] = -1;
    QuadrilateralCoo[2] =  3;
    QuadrilateralCoo[3] = -1;
    QuadrilateralCoo[4] =  3;
    QuadrilateralCoo[5] =  2;
    QuadrilateralCoo[6] = -1;
    QuadrilateralCoo[7] =  2;

    PointCoo[0] = 2.;
    PointCoo[1] = 0;

    printf("Starting to calculate the parametric co-ordinates\n");

    // calculate the parametric coordinates for this point
    F_InBilMap_Quadrilateral(QuadrilateralCoo, PointCoo, &xi1, &xi2, &both_roots_wanted);

    printf("Parametric Coordinates for this point are: xi1 = %.6f, xi2 = %.6f\n", xi1, xi2);

    double interpol_function_array[4]; // to store the interpolation function at this point
    F_Interpol_Function_Calculator_Quad(QuadrilateralCoo, PointCoo, interpol_function_array);
    printf("Interpolation function values at this point are:-\n");
    for (size_t i = 1; i <= 4; i++)
    {
        printf("Interpolation Function(%d) = %.6f\n", i, interpol_function_array[i - 1]);
    }

    // interpolate this quadrilateral using evaluated parametric co-ordinates
    double interpolated_x, interpolated_y;
    F_interpolate_quad(QuadrilateralCoo, xi1, xi2, &interpolated_x, &interpolated_y);
    printf("Interpolated point is (%f,%f)\n", interpolated_x, interpolated_y);
        

    /********************************************/
    /******* For both root ***********/
    /********************************************/
    both_roots_wanted = true;
    double both_xi1[2], both_xi2[2];

    QuadrilateralCoo[0] = -0.5;
    QuadrilateralCoo[1] = -0.5;
    QuadrilateralCoo[2] = 3;
    QuadrilateralCoo[3] = -1;
    QuadrilateralCoo[4] = 3;
    QuadrilateralCoo[5] = 2;
    QuadrilateralCoo[6] = -1;
    QuadrilateralCoo[7] = 2;

    PointCoo[0] = 2.;
    PointCoo[1] = 0;

    printf("................\n");
    printf("Starting to calculate the parametric co-ordinates\n");

    // calculate the parametric coordinates for this point
    F_InBilMap_Quadrilateral(QuadrilateralCoo, PointCoo, both_xi1, both_xi2, &both_roots_wanted);

    printf("  One Parametric Coordinates for this point are: xi1 = %.6f, xi2 = %.6f\n", both_xi1[0], both_xi2[0]);
    printf("Other Parametric Coordinates for this point are: xi1 = %.6f, xi2 = %.6f\n", both_xi1[1], both_xi2[1]);




}

