#include "materials/Nucleus.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::stringstream;
using std::setw;
using std::setfill;

#ifdef VECGEOM_NVCC
#ifndef VECGEOM_NVCC_DEVICE
#define fNuclei fNucleiHost
#define fIsoList fIsoListHost
#define fNatIsoList fNatIsoListHost
using std::find;
#else
template<class InputIterator, class T>
VECGEOM_CUDA_HEADER_BOTH
  InputIterator find (InputIterator first, InputIterator last, const T& val)
{
  while (first!=last) {
    if (*first==val) return first;
    ++first;
  }
  return last;
}
#define fNuclei fNucleiDev
#define fIsoList fIsoListDev
#define fNatIsoList fNatIsoListDev
#endif
#else
using std::find;
#include <string.h>
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECGEOM_NVCC
//________________________________________________________________________________________________
// Auxiliary functions to write the file

static void SwitchFile(std::stringstream &outline, unsigned int nfiles, unsigned int nfunc = 0)
{
  if (nfunc) {
    outline << "} // anonymous namespace\n\n";
    outline << "VECGEOM_CUDA_HEADER_BOTH" << endl << "void CreateParticle" << setfill('0') << setw(4) << nfiles-1 << "() {" << endl
            << setfill(' ') << setw(0);
    for (unsigned int i = 0; i < nfunc; ++i)
      outline << "  InternalCreateNuclei" << setfill('0') << setw(4) << i << "();" << endl;
  }
  outline << "}" << endl << endl;
  outline << "} // End of inline namespace" << endl;
  outline << "} // End of vecgeom namespace" << endl;
  outline << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
  outline << "#pragma clang optimize on" << endl;
  outline << "#endif" << endl;

  std::stringstream filename;
  filename << "CreateNuclei" << setfill('0') << setw(4) << nfiles-1 << ".cpp"
           << setfill(' ') << setw(0);

  std::ofstream outfile;
  outfile.open(filename.str());
  outfile << outline.str();
  outfile.close();

  outline.str("");
  outline.clear(); // Clear state flags.
}

static void StartFile(std::stringstream &outline, unsigned int nfiles, bool oneFuncPer = false)
{
  outline << "// This files was autogenerated by vecgeom::Nuclei::ReadFile\n\n";
  outline << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
  outline << "#pragma clang optimize off" << endl;
  outline << "#endif" << endl;
  outline << "#include \"materials/Nucleus.h\"" << endl;
  /*
  outline << "#ifdef VECGEOM_NVCC" << endl << "#include \"base/Vector.h\"" << endl
          << "template <typename T>" << endl << "using vector = vecgeom::Vector<T>;" << endl
          << "#else" << endl << "using std::vector;" << endl << "#endif" << endl << endl;
  */
  outline << "namespace vecgeom {" << endl;
  outline << "inline namespace VECGEOM_IMPL_NAMESPACE {" << endl << endl;
  outline << endl << "//" << setw(80) << setfill('_') << "_" << endl << setfill(' ') << setw(0);
  if (!oneFuncPer) {
    outline << "VECGEOM_CUDA_HEADER_BOTH" << endl << "void CreateNuclei" << setfill('0') << setw(4) << nfiles << "() {" << endl
            << setfill(' ') << setw(0);
    outline << "   Nucleus *nuc = 0;\n";
  } else {
    outline << "namespace {\n";
  }
}

#else
VECGEOM_CUDA_HEADER_BOTH
char *strncpy(char *dest, const char *src, size_t n)
{
  char *ret = dest;
  do {
    if (!n--) return ret;
  } while (*dest++ = *src++);
  while (n--)
    *dest++ = 0;
  return ret;
}
#endif

//________________________________________________________________________________________________

#ifdef VECGEOM_NVCC
class Nucleus;
VECGEOM_CUDA_HEADER_DEVICE vecgeom::map<int, Nucleus*> *fNucleiDev = nullptr; 
vecgeom::map<int, Nucleus*> *fNucleiHost = nullptr;
VECGEOM_CUDA_HEADER_DEVICE vecgeom::map<int, vector<Nucleus* >> *fIsoListDev = nullptr;
vecgeom::map<int, vector<Nucleus* >> *fIsoListHost = nullptr;
VECGEOM_CUDA_HEADER_DEVICE vecgeom::map<int, vector<Nucleus* >> *fNatIsoListDev = nullptr;
vecgeom::map<int, vector<Nucleus* >> *fNatIsoListHost = nullptr;
#else
std::map<int, Nucleus *> *Nucleus::fNuclei = nullptr;
std::map<int, std::vector<Nucleus *>> *Nucleus::fIsoList = nullptr;
std::map<int, std::vector<Nucleus *>> *Nucleus::fNatIsoList = nullptr;

std::ostream &operator<<(std::ostream &os, const Nucleus &nuc)
{
  os << nuc.Name() << ": A=" << nuc.fA;
  if (nuc.fNatab > 0) os << " Ab=" << nuc.fNatab;
  if (nuc.fIso != 0) os << " Level=" << nuc.fIsolevel << "MeV";
  if (nuc.fLife > 0) os << " Lifetime=" << nuc.fLife << "s";
  return os;
}
#endif

//________________________________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
Nucleus::Nucleus(const char *name, int n, int z, int iso, double a, double dm, double life, double natab, double toxa,
                 double toxb, int ind1, int ind2)
    : fN(n), fZ(z), fIso(iso), fA(a), fIsolevel(dm), fLife(life), fNatab(natab), fToxa(toxa), fToxb(toxb),
      fInd1(ind1), fInd2(ind2)
{

  strncpy(fName,name,49);
  fName[49]='\0';

  int zniso = 10000 * fZ + 10 * fN + fIso;

  if(!fNuclei) fNuclei = new map<int, Nucleus *>;
  if (fNuclei->count(zniso) != 0) {
     printf("Nucleus %d already there\n",zniso);
    return;
  }

  (*fNuclei)[zniso] = this;

  if(!fIsoList) fIsoList = new map<int, vector<Nucleus *>>;
  (*fIsoList)[fZ].push_back(this);
  if (natab > 0) {
    if(!fNatIsoList) fNatIsoList = new map<int, vector<Nucleus *>>;
    (*fNatIsoList)[fZ].push_back(this);
  }
}

#ifndef VECGEOM_NVCC
//________________________________________________________________________________________________
void Nucleus::ReadFile(std::string infilename, bool output)
{
  std::string name;
  int n, z, iso;
  double a, dm, life;
  int da, dz, diso;
  double br, qval, natab, toxa, toxb;
  int ind1, ind2;
  constexpr bool kOneFuncPer = false;
  constexpr unsigned int kMaxLines = 2400;

  // 100 real	1m11.904s user	1m6.876s  sys	0m4.693s
  // 500 real	0m20.706s user	0m19.036s sys	0m1.286s
  // 1000 real	0m10.877s user	0m10.190s sys	0m0.648s
  // 1500 real	0m8.175s user	0m7.669s sys	0m0.475s
  // 2000 real	0m7.349s user	0m6.895s sys	0m0.424s
  // 2400 real	0m6.286s user	0m5.904s sys	0m0.356s
  // 10000 real	0m4.664s user	0m4.418s sys	0m0.237s

  std::ofstream outfile;
  ifstream infile(infilename);
  std::string line;
  std::stringstream outline;
  std::stringstream filename;

  unsigned int nlines = 0;
  unsigned int nfiles = 0;
  unsigned int nfuncs = 0;
  unsigned int nnuclei= 0;

  getline(infile, line); // Get title
  Nucleus *nuc = 0;
  while (getline(infile, line)) {
    Getmat(line, n, z, iso, name, a, dm, life, da, dz, diso, br, qval, natab, toxa, toxb, ind1, ind2);
    if (z == 0) continue; // skip neutron

    int zniso = 10000 * z + 10 * n + iso;
    if (Nucleus::Nuclei().count(zniso) == 0) {
       nuc = new Nucleus(name.c_str(), n, z, iso, a, dm, life, natab, toxa, toxb, ind1, ind2);
      if (output) {
	 // Approximate number of lines 
	 if(nlines == 0 || nlines > kMaxLines) {
	    if(nnuclei > 0) {
	       SwitchFile(outline, nfiles, nfuncs);
	    }

	    StartFile(outline, nfiles, kOneFuncPer);
	    ++nfiles;

	    nlines = 0;
	    nfuncs = 0;
	 }
	 ++nnuclei;
        outline<< endl << "   // Adding " << nuc->Name() << endl;
        outline << "   nuc = new Nucleus("
                << "\"" << name << "\""
                << "," << n << "," << z << "," << iso << "," << a << "," << dm << "," << life << "," << natab << ","
                << toxa << "," << toxb << "," << ind1 << "," << ind2 << ");" << endl;
        ++nlines;
      }
    }
    if (da != 0 || dz != 0 || diso != 0) {
      nuc->AddDecay(da, dz, diso, qval, br);
      if (output) {
        outline << "   nuc->AddDecay(" << da << "," << dz << "," << diso << "," << qval << "," << br << ");" << endl;
	++nlines;
      }
    }
  }

  if (nlines != 0) {
     SwitchFile(outline, nfiles, nfuncs);
  }

  for (auto inuc = Nucleus::Nuclei().begin(); inuc != Nucleus::Nuclei().end(); ++inuc)
    inuc->second->NormDecay();

  if (output) {

    outfile.open("CreateNuclei.cpp");
    outfile << "// This files was autogenerated by vecgeom::Nuclei::ReadFile\n\n";
    outfile << "#include \"materials/Nucleus.h\"" << endl;
    outfile << "namespace vecgeom {" << endl;
    outfile << "inline namespace VECGEOM_IMPL_NAMESPACE {" << endl << endl;
    for (unsigned int i = 0; i < nfiles; ++i) {
      outfile << "VECGEOM_CUDA_HEADER_BOTH\n";
      outfile << "void CreateNuclei" << setfill('0') << setw(4) << i << "();" << endl;
     }
     
     outfile << "\n#ifdef VECGEOM_NVCC\n";
     outfile << "VECGEOM_CUDA_HEADER_DEVICE bool fgCreateNucleiInitDoneDev = false;\n";
     outfile << "#endif\n";
     
     outfile << endl << "//" << setw(80) << setfill('_') << "_" << endl << setfill(' ') << setw(0);
     outfile << "VECGEOM_CUDA_HEADER_BOTH" << endl << "void Nucleus::CreateNuclei() {" << endl;
     outfile << "#ifndef VECGEOM_NVCC\n";
     outfile << "  static bool fgCreateNucleiInitDone = false;\n";
     outfile << "#else\n";
     outfile << "  bool &fgCreateNucleiInitDone(fgCreateNucleiInitDoneDev);\n";
     outfile << "#endif\n";
     outfile << "  if (fgCreateNucleiInitDone) return;\n";
     outfile << "  fgCreateNucleiInitDone = true;\n";
     
     for (unsigned int i = 0; i < nfiles; ++i)
	outfile << "  CreateNuclei" << setfill('0') << setw(4) << i << "();" << endl;
     outfile << "} // End of CreateNuclei" << endl;
     outfile << "} // End of inline namespace" << endl;
     outfile << "} // End of vecgeom namespace" << endl;
     outfile << "#if defined(__clang__) && !defined(__APPLE__)" << endl;
     outfile << "#pragma clang optimize on" << endl;
     outfile << "#endif" << endl;
    outfile.close();

  }
}

void Nucleus::Getmat(std::string line, int &n, int &z, int &iso, std::string &name, double &a, double &dm, double &life,
                     int &da, int &dz, int &diso, double &br, double &qval, double &natab, double &toxa, double &toxb,
                     int &ind1, int &ind2)
{
  int beg     = 5;
  int ic      = 0;
  int len[17] = {5, 5, 5, 5, 15, 15, 15, 5, 5, 5, 15, 15, 15, 15, 15, 5, 5};

  beg = 5;
  ic  = 0;
  stringstream(line.substr(beg, len[ic])) >> n;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> z;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> iso;

  beg += len[ic];
  ++ic;
  name = line.substr(beg, len[ic]);
  name = name.substr(name.find_first_not_of(" "), name.find_last_not_of(" ") - name.find_first_not_of(" ") + 1);

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> a;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> dm;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> life;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> da;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> dz;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> diso;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> br;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> qval;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> natab;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> toxa;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> toxb;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> ind1;

  beg += len[ic];
  ++ic;
  stringstream(line.substr(beg, len[ic])) >> ind2;
}
#endif

//________________________________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
void Nucleus::NormDecay()
{
  double brt = 0;
  for (auto idec = fDecayList.begin(); idec != fDecayList.end(); ++idec)
    brt += idec->Br();
  brt = 1 / brt;
  for (auto idec = fDecayList.begin(); idec != fDecayList.end(); ++idec)
    idec->Br(100 * idec->Br() * brt);
}

//________________________________________________________________________________________________
VECGEOM_CUDA_HEADER_BOTH
void Nucleus::AddDecay(int da, int dz, int diso, double qval, double br)
{
  Decay dec(da, dz, diso, qval, br);
  bool found = false;
  if (find(fDecayList.begin(), fDecayList.end(), dec) != fDecayList.end()) {
    printf("Decay already there!\n");
    found = true;
  }
  if (!found) fDecayList.push_back(dec);
}


//________________________________________________________________________________________________
const std::string Nucleus::Decay::Name() const
{
  std::stringstream name;
  name << "(" << fBr << "%) ";
  if (fDz == -2 && fDa == -4) {
    name << "Alpha";
    if (fDiso != 0) name << " iso";
  } else if (fDz == 1 && fDa == 0) {
    name << "Beta-";
    if (fDiso != 0) name << " iso";
  } else if (fDz == -1 && fDa == 0) {
    name << "Beta+";
    if (fDiso != 0) name << " iso";
  } else if (fDz == 0 && fDa == 0 && fDiso == -1)
    name << "IC";
  else if (fDz == 1000)
    name << "Fission";
  else
    name << fDa << ":" << fDz << ":" << fDiso;
  return name.str();
}
}
}
