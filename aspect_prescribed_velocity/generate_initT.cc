#include <iostream>
#include <math.h>


int main()
{

    const double diffusivity = 1.0e-6;
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
    const double xt = 100.0e3;
    const double x0 = 150.0e3;
            
    const int nx = 701;
    const int nz = 351;
    const double xmin = 0.0e3;
    const double xmax = 700.0e3;
    const double zmin = 0.0e3;
    const double zmax = 350.0e3;
  
    for( int j = 0; j < ny; ++j ) {
        for (int i = 0; i < nx; ++i) {
            double xn = i*(xmax-xmin)/(nx-1);
            double zn = j*(zmax-zmin)/(nz-1);
            double yn = 0.0;
            double age = 0.0;
            double depth  = zmax - zn;
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
