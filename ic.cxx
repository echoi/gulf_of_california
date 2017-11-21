#include <iostream>
#include <math.h>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"

#include "ic-read-temp.hpp"
#include "ic.hpp"

namespace {

    class Zone
    {
    public:
        virtual ~Zone() {};
        virtual bool contains(const double x[NDIMS]) const = 0;
    };

    class Empty_zone : public Zone
    {
    public:
        bool contains(const double x[NDIMS]) const {return false;}
    };


    class Planar_zone : public Zone
    {
    private:
        const double az, incl;
        const double halfwidth; // in meter
#ifdef THREED
        const double ymin, ymax; // in meter
#endif
        const double zmin, zmax; // in meter
        const double *x0;

    public:
        Planar_zone(const double center[NDIMS], double azimuth, double inclination, double halfwidth_,
#ifdef THREED
                    double ymin_, double ymax_,
#endif
                    double zmin_, double zmax_) :
            az(std::tan(azimuth * DEG2RAD)), incl(1/std::tan(inclination * DEG2RAD)), halfwidth(halfwidth_),
#ifdef THREED
            ymin(ymin_), ymax(ymax_),
#endif
            zmin(zmin_), zmax(zmax_),
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {}

        bool contains(const double x[NDIMS]) const
        {
            // Is x within halfwidth distance to a plane cutting through x0?
            return (x[NDIMS-1] > zmin &&
                    x[NDIMS-1] < zmax &&
#ifdef THREED
                    x[1] > ymin &&
                    x[1] < ymax &&
#endif
                    std::fabs( (x[0] - x0[0])
#ifdef THREED
                               - az * (x[1] - x0[1])
#endif
                               + incl * (x[NDIMS-1] - x0[NDIMS-1]) ) < halfwidth );
        }
    };


    class Ellipsoidal_zone : public Zone
    {
    private:
        const double *x0;
        double semi_axis2[NDIMS];

    public:
        Ellipsoidal_zone(const double center[NDIMS], const double semi_axis[NDIMS]) :
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {
            for(int i=0; i<NDIMS; i++)
                semi_axis2[i] =  semi_axis[i] * semi_axis[i];
        }

        bool contains(const double x[NDIMS]) const
        {
            return ( (x[0] - x0[0])*(x[0] - x0[0])/semi_axis2[0]
#ifdef THREED
                     + (x[1] - x0[1])*(x[1] - x0[1])/semi_axis2[1]
#endif
                     + (x[NDIMS-1] - x0[NDIMS-1])*(x[NDIMS-1] - x0[NDIMS-1])/semi_axis2[NDIMS-1] < 1 );
        }
    };

} // anonymous namespace


void initial_stress_state(const Param &param, const Variables &var,
                          tensor_t &stress, double_vec &stressyy, tensor_t &strain,
                          double &compensation_pressure)
{
    if (param.control.gravity == 0) {
        compensation_pressure = 0;
        return;
    }

    // lithostatic condition for stress and strain
    double rho = var.mat->rho(0);
    double ks = var.mat->bulkm(0);

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        double p = ref_pressure(param, zcenter);
        if (param.control.ref_pressure_option == 1 ||
            param.control.ref_pressure_option == 2) {
            ks = var.mat->bulkm(e);
        }

        for (int i=0; i<NDIMS; ++i) {
            stress[e][i] = -p;
            strain[e][i] = -p / ks / NDIMS;
        }
        if (param.mat.is_plane_strain)
            stressyy[e] = -p;
    }

    compensation_pressure = ref_pressure(param, -param.mesh.zlength);
}


void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain)
{
    Zone *weakzone;

    // TODO: adding different types of weak zone
    double plane_center[NDIMS]; // this variable must outlive weakzone
    switch (param.ic.weakzone_option) {
    case 0:
        weakzone = new Empty_zone();
        break;
    case 1:
        // a planar weak zone, cut through top center
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        weakzone = new Planar_zone(plane_center,
                                   param.ic.weakzone_azimuth,
                                   param.ic.weakzone_inclination,
                                   param.ic.weakzone_halfwidth * param.mesh.resolution,
#ifdef THREED
                                   param.ic.weakzone_y_min * param.mesh.ylength,
                                   param.ic.weakzone_y_max * param.mesh.ylength,
#endif
                                   -param.ic.weakzone_depth_max * param.mesh.zlength,
                                   -param.ic.weakzone_depth_min * param.mesh.zlength);
        break;
    case 2:
        // a ellipsoidal weak zone
        double semi_axis[NDIMS];
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
        semi_axis[0] = param.ic.weakzone_xsemi_axis;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
        semi_axis[1] = param.ic.weakzone_ysemi_axis;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        semi_axis[NDIMS-1] = param.ic.weakzone_zsemi_axis;
        weakzone = new Ellipsoidal_zone(plane_center, semi_axis);
        break;
    default:
        std::cerr << "Error: unknown weakzone_option: " << param.ic.weakzone_option << '\n';
        std::exit(1);
    }

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        // the coordinate of the center of this element
        double center[NDIMS] = {0};
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int d=0; d<NDIMS; ++d) {
                center[d] += (*var.coord)[conn[i]][d];
            }
        }
        for (int d=0; d<NDIMS; ++d) {
            center[d] /= NODES_PER_ELEM;
        }

        if (weakzone->contains(center))
            plstrain[e] = param.ic.weakzone_plstrain;

        // Find the most abundant marker mattype in this element
        // int_vec &a = (*var.elemmarkers)[e];
        // int material = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
    }

    delete weakzone;
}


void initial_temperature(const Param &param, const Variables &var,
                         double_vec &temperature)
{
    switch(param.ic.temperature_option) {
    case 0:
        {
            const double age = param.ic.oceanic_plate_age_in_yr * YEAR2SEC;
            const MatProps &mat = *var.mat;
            const double diffusivity = mat.k(0) / mat.rho(0) / mat.cp(0); // thermal diffusivity of 0th element

            for (int i=0; i<var.nnode; ++i) {
                double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * age);
                temperature[i] = param.bc.surface_temperature +
                    (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
            }
            break;
        }
    case 1:
        {
            const MatProps &mat = *var.mat;
            const double diffusivity = mat.k(0) / mat.rho(0) / mat.cp(0); // thermal diffusivity of 0th element

            const double dip_FAR = 15.0 * M_PI / 180.0;
            const double tan_delta = tan(dip_FAR);
            const double tan_delta_sqrd = tan_delta * tan_delta;
            const double strike_Ridge = 7.0 * M_PI / 180.0;
            const double tan_theta = tan(strike_Ridge);
            const double cos_del_cos_theta = cos( dip_FAR * cos(strike_Ridge) );
            const double obliquity_FAR = 20.0 * M_PI / 180.0;
            const double cos_alpha_FAR = cos(obliquity_FAR);
            const double obliquity_PAC = 330.0 * M_PI / 180.0;
            const double cos_alpha_PAC = cos(obliquity_PAC);
            const double v_half = 0.106 / YEAR2SEC; // m/yr to m/s.
            const double xt = 200.0e3;
            const double x0 = 150.0e3;
            
            for (int i=0; i<var.nnode; ++i) {
                double xn = (*var.coord)[i][0];
                double yn = (*var.coord)[i][1];
                double zn = (*var.coord)[i][2];
                double age = 0.0;
                double depth  = std::fabs(zn);
                // If in or beneath unsubducted plates
                if( xn <= xt ) {
                    double xr = (x0 + tan_theta * (tan_theta*xn-yn))/(1.0+tan_theta*tan_theta);
                    double yr = tan_theta * (xr-xn) + yn;
                    double dist = sqrt( (xr-xn)*(xr-xn) + (yr-yn)*(yr-yn) );
                    age = dist / (v_half * cos_alpha);
                }
                // If in the NA plate
                else if( xn > xt && zn > -(xn-xt) * tan_delta ) {
                    age = 100.0e6 * YEAR2SEC; // 100 My in seconds
                }
                // If in or beneath the subducted plate
                else if( xn > xt && zn <= -(xn-xt) * tan_delta ) {
                    double xs = (xn - tan_delta * zn + tan_delta_sqrd * xt) / (1.0 + tan_delta_sqrd );
                    double ys = yn;
                    double zs = -tan_delta * ( xs - xt );
                    double xr = (x0 + tan_theta * (tan_theta*xs-ys))/(1.0+tan_theta*tan_theta);
                    double yr = tan_theta * (xr-xs) + ys;
                    double yt = tan_theta * (xt-xs) + ys;
                    d_tr = sqrt( (xr-xt)*(xr-xt) + (yr-yt)*(yr-yt) );
                    d_st = sqrt( (xt-xs)*(xt-xs) + (yt-ys)*(yt-ys) );
                    age = d_tr / (v_half * cos_alpha) + d_st / ( v_half * cos_alpha * cos_del_cos_theta);
                    depth = sqrt( (xs-xn)*(xs-xn) + (zs-zn)*(zs-zn) );
                }
            
                double w = 0.5 * depth / std::sqrt( diffusivity * age );
                temperature[i] = param.bc.surface_temperature +
                    (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
            }
            break;
        }

    case 90:
        read_external_temperature_from_comsol(param, var, *var.temperature);
        break;
    default:
        std::cout << "Error: unknown ic.temperature option: " << param.ic.temperature_option << '\n';
        std::exit(1);
    }
}


