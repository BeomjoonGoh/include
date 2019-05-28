#ifndef CONSTANTS_H
#define CONSTANTS_H

// Values given are in SI unit.  Use the following constants by redefining them with shortened names within the code.
// For example,
//  double c_0 = Constants::speed_of_light_in_vacuum;

namespace Constants
{
  //                | Name                               | Value               | Unit
  //----------------+------------------------------------+---------------------+---------------

  // Universal Constants 
  extern const double speed_of_light_in_vacuum            {2.99792458e+8   };   // m.s-1       
  extern const double reduced_planck_constant             {1.054571726e-34 };   // J.s         
  extern const double planck_constant                     {6.62606957e-34  };   // J.s         
  extern const double newtonian_constant_of_gravitation   {6.67384e-11     };   // m3.kg-1.s-2 

  // Natural Units 
  extern const double planck_length                       {1.616199e-35    };   // m           
  extern const double planck_mass                         {2.17651e-8      };   // kg          
  extern const double planck_time                         {5.39106e-44     };   // s           
  extern const double planck_charge                       {1.875545956e-18 };   // C           
  extern const double planck_temperature                  {1.416833e+32    };   // K           

  // Electromagnetic Constants 
  extern const double vacuum_permeability                 {1.256637061e-6  };   // N.A-2       
  extern const double vacuum_permittivity                 {8.854187817e-12 };   // F.m-1       
  extern const double elementary_charge                   {1.602176565e-19 };   // C           
  extern const double bohr_magneton                       {9.27400968e-24  };   // J.T-1       
  extern const double nuclear_magneton                    {5.05078353e-27  };   // J.T-1       
  extern const double conductance_quantum                 {7.7480917346e-5 };   // S           
  extern const double magnetic_conductance_quantum        {2.067833758e-15 };   // Wb          

  // Atomic & Nuclear Constants 
  extern const double bohr_radius                         {5.2917721092e-11};   // m           
  extern const double electron_mass                       {9.10938291e-31  };   // kg          
  extern const double fine_structure_constant             {7.2973525698e-3 };
  extern const double proton_mass                         {1.672621777e-27 };   // kg          
  extern const double rydberg_constant                    {1.09737315685e+7};   // m-1         

  // Physico-Chemical Constants 
  extern const double atomic_mass_constant                {1.660538921e-27 };   // kg          
  extern const double avogadro_number                     {6.02214129e+23  };   // mol-1       
  extern const double boltzmann_constant                  {1.3806488e-23   };   // J.K-1       
  extern const double faraday_constant                    {9.64853365e+4   };   // C.mol-1     
  extern const double gas_constant                        {8.3144621       };   // J.K-1.mol-1 
  extern const double stefan_boltzmann_constant           {5.670373e-8     };   // W.m-2.K-4   

  // Adopted values 
  extern const double standard_acceleration_of_gravity    {9.80665         };   // m.s-2       
  extern const double standard_atmosphere                 {1.01325e+5      };   // Pa          

}
#endif
