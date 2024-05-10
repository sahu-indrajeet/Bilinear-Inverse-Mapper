# [[Preprint]](https://dx.doi.org/10.2139/ssrn.4790071)

# Inverse Mapping Algorithm for Bilinear Quadrilateral
A powerful C library for finding the parametric co-ordinates of bilinear mapping for any given point with respect to a quadrilateral. These parametric co-ordinates are also utilised in calculating the interpolation functions or shape functions at the given point. The library also helps in evaluating the Barycentric co-ordinates of a point in case of a triangle.


Given a point and a quadrilateral in the physical space, the function $\phi^{-1}(x)$ defines the inverse mapping to a corresponding point in the parametric space. The algorithm calculates the parametric co-ordinates $(\xi_1,\xi_2)$ corresponding to the given point $(x,y)$ in the physical space through this inverse mapping.

![Mapping between parametric co-ordinates and spatial coordinates](img/Quad_Parametric_to_Physical_Mapping.svg)

## Installation
To install the library, simply clone this repository and include the header files in your project. Make sure to link the library during compilation.

Refer to examples in this repository for a quick guide on this library's utility. 
