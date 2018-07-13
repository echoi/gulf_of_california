#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>

int main()
{

    const double diffusivity = 1.0e-6;
    const double dip_FAR = 15.0 * M_PI / 180.0;
    const double tan_delta = tan(dip_FAR);
    const double cos_delta = cos(dip_FAR);
    const double tan_delta_sqrd = tan_delta * tan_delta;
    const double strike_Ridge = 7.0 * M_PI / 180.0;
    const double tan_theta = tan(strike_Ridge);
    const double cos_del_cos_theta = cos( dip_FAR * cos(strike_Ridge) );
    const double obliquity_FAR = 20.0 * M_PI / 180.0;
    const double cos_alpha_FAR = cos(obliquity_FAR);
    const double obliquity_PAC = 330.0 * M_PI / 180.0;
    const double cos_alpha_PAC = cos(obliquity_PAC);
    const double YEAR2SEC = 3.1536e7;
    const double v_half = 0.03 / YEAR2SEC; // m/yr to m/s.

    const double xt = 100.0e3;
    const double x0 = 0.0e3;
    const double surface_temperature = 293.0;
    const double mantle_temperature = 1593.0;
            
    const int nx = 701;
    const int nz = 351;
    const double xmin = 0.0e3;
    const double xmax = 700.0e3;
    const double zmin = 0.0e3;
    const double zmax = 350.0e3;
    const int num_phase = 4;
  
    std::ofstream initC_file;
    initC_file.open("goc_initC_Baja.txt");
    initC_file << "# POINTS: " << nx <<" "<< nz << std::endl;
    initC_file << "# Columns: x y NA_upper_crust NA_lower_crust NA_mantle" << std::endl;
    for( int j = 0; j < nz; ++j ) {
        for (int i = 0; i < nx; ++i) {
            std::vector<double> compositions(num_phase,0.0);
            double xn = i*(xmax-xmin)/(nx-1);
            double zn = j*(zmax-zmin)/(nz-1);
            // If in or beneath unsubducted plates, all the comp fields are zero. So don't do anything.
            // If in the NA plate
            if( xn > xt && (zmax - zn) < ( xn - xt ) * tan_delta ) {
                // if in the crust
                if(zn > (zmax - 35.0e3)) {
                    // if in the mafic Baja crust
                    if((xn - xt) < 100.e3) {
                        compositions[0] = 1.0;
                        std::cerr << "(0) " << xn <<" "<< zn <<" "<< compositions[0]<<" "<< compositions[1] <<" "<< compositions[2] <<" "<< compositions[3] << std::endl;
                    }
                    // if in the felsec Baja crust
                    else if((xn - xt) < 200.e3) {
                        compositions[1] = 1.0;
                        std::cerr << "(1) " << xn <<" "<< zn <<" "<< compositions[0]<<" "<< compositions[1] <<" "<< compositions[2] <<" "<< compositions[3] << std::endl;
                    }
                    // if in the NA crust
                    else {
                        compositions[2] = 1.0;
                        std::cerr << "(2) " << xn <<" "<< zn <<" "<< compositions[0]<<" "<< compositions[1] <<" "<< compositions[2] <<" "<< compositions[3] << std::endl;
                    }
                }
                // if in the NA mantle
                else {
                    compositions[3] = 1.0;
                    std::cerr << "(3) " << xn <<" "<< zn <<" "<< compositions[0]<<" "<< compositions[1] <<" "<< compositions[2] <<" "<< compositions[3] << std::endl;
                }
            }
            initC_file << xn <<" "<< zn;
            for (std::vector<double>::iterator it = compositions.begin(); it != compositions.end(); ++it)
                initC_file << ' ' << *it;
            initC_file << std::endl;
            //std::cerr << xn <<" "<< zn <<" "<< compositions[0]<<" "<< compositions[1] <<" "<< compositions[2] << std::endl;
        }
    }
    initC_file.close();

    return 0;
}
