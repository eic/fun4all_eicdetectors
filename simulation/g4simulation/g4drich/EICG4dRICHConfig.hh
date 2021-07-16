#ifndef EICG4dRICHConfig_HH
#define EICG4dRICHConfig_HH

#include <G4GenericMessenger.hh>
#include <G4SystemOfUnits.hh>

struct EICG4dRICHConfig 
{

  // define here Global volume parameters
  double RIn = 9.;
  double ROut = 220 * cm;
  double ThicknessZ = 170 * cm;
  double PosZ = 200. * cm;

  //
  // RICH is made of 6 identical trucated cone sectors (TBD: probably one or
  // more sectors shall be different due to the asymmetry in the beam pipes)
  // Each sector components have identical relative disposition
  // Each sector is symmetric respect to the cone bysector ...

  /* model data:
          material, geometry and placement are defined here
   */

  G4String model_file = "drich-model0.txt";

  /*
    double vessel_z0 = 300*cm;  // position of the vessel entrance window along
    z respect to the overall spectrometer ??? double vessel_dz = 170*cm; //
    Vessel longitudinal length (external) - 0 thickness for the moment double
    vessel_radius_in = 10*cm; // Vessel internal radius (around beam pipe)
    double vessel_radius_out0 = 130.*cm; // Vessel external radius at entrance
    double vessel_radius_out1 = 210.*cm; // Vessel external radius at exit
    (behind the mirror) double aerogel_dz = 2*cm; // aerogel thickness double
    filter_dz = 0.3*mm; // wavelength filter thickness double mirror_radius =
    290*cm; // mirror radius double mirror_center_z = -210*cm; // mirror center
    z position - relative to vessel center double mirror_center_t = 145*cm; //
    mirror center transverse position - relative to vessel center double
    mirror_theta = 40*deg; // mirror polar aperture

    double sensor_z = 0; // sensor center along z - relative to vessel center
    double sensor_t = 0; // sensor center along tranverse z - relative to vessel
    center
  */

  /* optical data:
         aerogel dispersione relation n_aerogel(lambda)
         aerogel scattering length vs lambda
         aerogel absorption length vs lambda
         wavelength filter transmission curve vs lambda
         gas dispersione relation n_gas(lambda)
         gas absorption lenghth vs lambda (likely constant)
         mirror refractive curve vs lambda
         sensor quantum efficiency (embedded in the code or external ?)
  */
  //  G4String optical_file = "drich-optical.txt";

  // optical properties
  int aerOptModel = 3; // aerogel optical model used to estimate the refractive Index
  double filter_thr = 300 * nm; // wavelength filter cutoff

  double aerogel_n = 1.03; // aerogel refractive index, at reference wavelength
  double gas_n = 1.0008;   // gas refractive index, at reference wavelength
  double wavelength_ref = 400 * nm; // reference wavelength (for aerogel, gas, mirror ...)

  // Create a global messenger that will be used
  inline EICG4dRICHConfig() 
  {
    static G4GenericMessenger *Messenger;

    if (!Messenger) 
    {
      // Set geant options
      Messenger = new G4GenericMessenger(this, "/EICG4dRICH/");

      //    Messenger->SetDirectory("Vessel"); // TBC if ok, we can define a dir
      //    for each sub-components of the EICG4dRICH
      /*
        auto cmd = Messenger->DeclarePropertyWithUnit("VesselLength", "cm",
        vessel_dz, "Rich Vessel Longitudinal size, along beam (z), cone
        height"); cmd.SetParameterName("vl",true); cmd.SetRange("vl>0,");
        cmd.SetDefaultValue("170");
        Messenger->DeclarePropertyWithUnit("VesselPositionZ0", "cm", vessel_z0,
        "Rich Vessel Entrance Window Position along z");
        Messenger->DeclarePropertyWithUnit("VesselRadiusIn", "cm",
        vessel_radius_in, "Rich Vessel Internal Radius (near beam pipe)");
        Messenger->DeclarePropertyWithUnit("VesselRadiusOut0", "cm",
        vessel_radius_out0, "Rich Vessel External Radius at entrance");
        Messenger->DeclarePropertyWithUnit("VesselRadiusOut1", "cm",
        vessel_radius_out1, "Rich Vessel External Radius at exit, after
        mirror"); Messenger->DeclarePropertyWithUnit("AerogelThick", "cm",
        aerogel_dz, "Aerogel Thickness");
        Messenger->DeclarePropertyWithUnit("FilterThick", "mm", filter_dz,
        "High-Pass Filter Thickness");
        Messenger->DeclarePropertyWithUnit("MirrorRadius", "cm", mirror_radius,
        "Mirror Radius"); Messenger->DeclarePropertyWithUnit("MirrorCenterZ",
        "cm", mirror_center_z, "Mirror Center - longitudinal, Z relative to
        Vessel center"); Messenger->DeclarePropertyWithUnit("MirrorCenterT",
        "cm", mirror_center_t, "Mirror Center (x/y plane) Transverse, relative
        to Vessel center"); Messenger->DeclarePropertyWithUnit("MirrorTheta",
        "deg", mirror_theta, "Mirror Polar Aperture");
        Messenger->DeclarePropertyWithUnit("SensorZ", "cm", sensor_z, "Sensor
        Center - longitudinal, Z relative to Vessel center");
        Messenger->DeclarePropertyWithUnit("SensorT", "cm", sensor_t, "Sensor
        Center - Transverse (x/y plane) relative to Vessel center");
      */

      Messenger->DeclareProperty("ModelFile", model_file,"Full path of the drich model text data");

      auto cmd0 = Messenger->DeclareProperty(
          "AerogelRefractiveIndex", aerogel_n,
          "Aerogel Radiator Refractive Index at reference wavelength");
      cmd0.SetParameterName("aerogeln", true);
      cmd0.SetRange("aerogeln>1.0");
      cmd0.SetDefaultValue("1.03");

      auto cmd1 = Messenger->DeclareProperty(
          "GasRefractiveIndex", gas_n,
          "Gas Radiator Refractive Index at reference wavelength");
      cmd1.SetParameterName("gasn", true);
      cmd1.SetRange("gasn>1.0");
      cmd1.SetDefaultValue("1.0008");

      auto cmd2 = Messenger->DeclarePropertyWithUnit(
          "ReferenceWavelength", "nm", wavelength_ref, "Reference wavelength");
      cmd2.SetParameterName("wlref", true);
      cmd2.SetRange("wlref>0");
      cmd2.SetDefaultValue("400.e-9");

      auto cmd3 = Messenger->DeclarePropertyWithUnit(
          "FilterThreshold", "nm", filter_thr,
          "Wavelength threshold of acrylic filter");
      cmd3.SetParameterName("wlfil", true);
      cmd3.SetRange("wlfil>0.");
      cmd3.SetDefaultValue("300.e-9");

      auto cmd4 = Messenger->DeclareProperty(
          "AerogelOpticalModel", aerOptModel,
          "Aerogel Optical Model for Refractive Index (0:Vorobiev, 1:Sellmeier "
          "CLAS12, 2: LHCb, 3: scale exp. points");
      cmd4.SetParameterName("aerOptMod", true);
      cmd4.SetDefaultValue(3);

      //      Messenger->DeclareProperty("GeometryFile", geometry_file, "Full
      //      path of Geometry configuration file (e.g. sensor array coordinates
      //      ...)");
      // Messenger->DeclareProperty("OpticalFile", optical_file, "Full path of
      // Optical configuration file (e.g. sensor quantum efficienct curve,
      // mirror transmittance, radiators dispersion relations ...)");
    }
  }
};

#endif // EICG4dRICHConfig_HH
