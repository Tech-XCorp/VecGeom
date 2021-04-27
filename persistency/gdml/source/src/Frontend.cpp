//!    \file Frontend.cpp
//!    \brief implements loading GDML to VecGeom
//!
//!    \authors Author:  Dmitry Savin <sd57@protonmail.ch>
//!

#include "Frontend.h"
#include "Backend.h"
#include "VecGeom/management/GeoManager.h"

#include <string>

namespace vgdml {

namespace Frontend {
std::unique_ptr<Middleware> Load(std::string const &aFilename, bool validate, double mm_unit, bool verbose)
{
  auto aBackend      = vgdml::Backend(validate);
  auto const aDOMDoc = aBackend.Load(aFilename);
  if (!aDOMDoc) {
    std::cerr << "== Error: GDML file " << aFilename << " could not be loaded\n";
    return nullptr;
  }
  vecgeom::GeoManager::SetMillimeterUnit((vecgeom::Precision)mm_unit);
  if (verbose == 1) std::cout << "(II) vgdml::Frontend::Load: VecGeom millimeter is " << mm_unit << "\n";
  auto aMiddleware            = new vgdml::Middleware();
  auto const loadedMiddleware = aMiddleware->Load(aDOMDoc);
  if (!loadedMiddleware) {
    delete aMiddleware;
    aMiddleware = nullptr;
  }
  return std::unique_ptr<Middleware>(aMiddleware);
}

} // namespace Frontend

} // namespace vgdml
