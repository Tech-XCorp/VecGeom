#include <VecCore/VecCore>
#include "base/Math.h"
#include "base/Vector3D.h"
#include "volumes/Tessellated.h"
#include "volumes/Extruded.h"

using namespace vecgeom;

vecgeom::UnplacedExtruded *ExtrudedMultiLayer()
{
  const size_t nvert             = 8;
  const size_t nsect             = 4;
  vecgeom::XtruVertex2 *vertices = new vecgeom::XtruVertex2[nvert];
  vecgeom::XtruSection *sections = new vecgeom::XtruSection[nsect];

  vertices[0].x = -3;
  vertices[0].y = -3;
  vertices[1].x = -3;
  vertices[1].y = 3;
  vertices[2].x = 3;
  vertices[2].y = 3;
  vertices[3].x = 3;
  vertices[3].y = -3;
  vertices[4].x = 1.5;
  vertices[4].y = -3;
  vertices[5].x = 1.5;
  vertices[5].y = 1.5;
  vertices[6].x = -1.5;
  vertices[6].y = 1.5;
  vertices[7].x = -1.5;
  vertices[7].y = -3;

  sections[0].fOrigin.Set(-2, 1, -4.0);
  sections[0].fScale = 1.5;
  sections[1].fOrigin.Set(0, 0, 1.0);
  sections[1].fScale = 0.5;
  sections[2].fOrigin.Set(0, 0, 1.5);
  sections[2].fScale = 0.7;
  sections[3].fOrigin.Set(2, 2, 4.0);
  sections[3].fScale = 0.9;

  UnplacedExtruded *xtru = new UnplacedExtruded(nvert, vertices, nsect, sections);
  return xtru;
}

size_t TessellatedOrb(double r, int ngrid, vecgeom::UnplacedTessellated &tsl)
{
  // Create a tessellated orb divided in ngrid*ngrid theta/phi cells
  // Sin/Cos tables
  std::vector<double> sth, cth, sph, cph;
  double dth = vecgeom::kPi / ngrid;
  double dph = vecgeom::kTwoPi / ngrid;

  for (int i = 0; i <= ngrid; ++i) {
    sth.push_back(vecCore::math::Sin(i * dth));
    cth.push_back(vecCore::math::Cos(i * dth));
    sph.push_back(vecCore::math::Sin(i * dph));
    cph.push_back(vecCore::math::Cos(i * dph));
  }

  auto Vtx = [&](int ith, int iph) {
    return vecgeom::Vector3D<double>(r * sth[ith] * cph[iph], r * sth[ith] * sph[iph], r * cth[ith]);
  };

  for (int ith = 0; ith < ngrid; ++ith) {
    for (int iph = 0; iph < ngrid; ++iph) {
      // First/last rows - > triangles
      if (ith == 0) {
        tsl.AddTriangularFacet(vecgeom::Vector3D<double>(0, 0, r), Vtx(ith + 1, iph), Vtx(ith + 1, iph + 1));
      } else if (ith == ngrid - 1) {
        tsl.AddTriangularFacet(Vtx(ith, iph), vecgeom::Vector3D<double>(0, 0, -r), Vtx(ith, iph + 1));
      } else {
        tsl.AddQuadrilateralFacet(Vtx(ith, iph), Vtx(ith + 1, iph), Vtx(ith + 1, iph + 1), Vtx(ith, iph + 1));
      }
    }
  }
  return (tsl.GetNFacets());
}
