#ifdef __CINT__

#pragma link C++ class std::pair<EICPIDDefs::PIDCandidate, EICPIDDefs::PIDDetector> +;
#pragma link C++ class std::pair<std::pair<EICPIDDefs::PIDCandidate, EICPIDDefs::PIDDetector>, float> +;
#pragma link C++ class std::map<std::pair<EICPIDDefs::PIDCandidate, EICPIDDefs::PIDDetector>, float> +;
#pragma link C++ class EICPIDParticlev1 +;

#endif /* __CINT__ */
