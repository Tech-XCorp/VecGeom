#ifndef VECGEOM_MANAGEMENT_CUDAMANAGER_H_
#define VECGEOM_MANAGEMENT_CUDAMANAGER_H_

#include <set>
#include <map>
#include "volumes/box.h"

namespace vecgeom {

class CudaManager {

private:

  bool synchronized;

  std::set<VUnplacedVolume const*> unplaced_volumes;
  std::set<LogicalVolume const*> logical_volumes;
  std::set<VPlacedVolume const*> placed_volumes;
  std::set<TransformationMatrix const*> matrices;
  std::set<Container<VPlacedVolume const*> *> daughters;

  typedef void const* CpuAddress;
  typedef void* GpuAddress;
  typedef std::map<CpuAddress, GpuAddress> MemoryMap;

  MemoryMap memory_map;

public:

  static CudaManager& Instance() {
    static CudaManager instance;
    return instance;
  }

  void LoadGeometry(LogicalVolume const *const volume);

  void Synchronize();

  void PrintContent() const;

private:

  CudaManager() {
    synchronized = true;
  }
  CudaManager(CudaManager const&);
  CudaManager& operator=(CudaManager const&);

  void CleanGpu();

  void ScanGeometry(LogicalVolume const *const volume);

  /**
   * Allocates all objects retrieved by ScanGeometry() on the GPU, storing
   * pointers in the memory table for future reference.
   */
  void AllocateGeometry();

  template <typename Type>
  CpuAddress ToCpuAddress(Type const *const ptr) {
    return static_cast<CpuAddress>(ptr);
  }

  template <typename Type>
  GpuAddress ToGpuAddress(Type *const ptr) {
    return static_cast<GpuAddress>(ptr);
  }

};

} // End namespace vecgeom

#endif // VECGEOM_MANAGEMENT_CUDAMANAGER_H_