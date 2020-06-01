/// \file NavIndexTable.cpp
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#include "VecGeom/management/NavIndexTable.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

NavIndex_t BuildNavIndexVisitor::apply(NavStatePath *state, int level, NavIndex_t mother, int dind)
{
  bool cacheTrans       = true;
  NavIndex_t new_mother = fCurrent;
  auto lv               = state->Top()->GetLogicalVolume();
  size_t nd             = lv->GetDaughters().size();
  // Check if matrix has to be cached for this node
  if (fLimitDepth > 0 && level > fLimitDepth && !lv->IsReqCaching()) cacheTrans = false;

  if (fValidate) {
    return NavIndexTable::Instance()->ValidateState(state);
  }
  // Size in bytes of the current node data
  size_t current_size = 3 * sizeof(unsigned int) + int(cacheTrans) * 12 * sizeof(double) + nd * sizeof(unsigned int);
  if (fDoCount) {
    fTableSize += current_size;
    return 0;
  }

  // Add data for the current element.

  // Fill the mother index for the current node
  fNavInd[fCurrent] = mother;

  // Fill the node index in the mother list of daughters
  if (mother > 0) fNavInd[mother + 3 + dind] = fCurrent;

  // Physical volume index
  fNavInd[fCurrent + 1] = (level >= 0) ? state->ValueAt(level) : 0;

  // Write current level in next byte
  auto content_ddt = (unsigned char *)(&fNavInd[fCurrent + 2]);
  assert(level < std::numeric_limits<unsigned char>::max() && "fatal: geometry deph more than 255 not supported");
  *content_ddt = (unsigned char)level;

  // Write number of daughters in next 2 bytes
  auto content_nd = (unsigned short *)(content_ddt + 1);
  assert(nd < std::numeric_limits<unsigned short>::max() && "fatal: not supporting more than 65535 daughters");
  *content_nd = (unsigned short)nd;

  // Write the flag if matrix is stored in the next byte
  auto content_hasm = (unsigned char *)(content_ddt + 3);
  *content_hasm     = 0;

  // Prepare the space for the daughter indices
  auto content_dind = &fNavInd[fCurrent + 3];
  for (size_t i = 0; i < nd; ++i)
    content_dind[i] = 0;

  fCurrent += 3 + nd;

  if (!cacheTrans) return new_mother;

  Transformation3D mat;
  // encode has_trans, translation and rotation flags in the content_hasm byte
  state->TopMatrix(mat);
  *content_hasm = 0x04 + 0x02 * (unsigned short)mat.HasTranslation() + (unsigned short)mat.HasRotation();

  // Write the transformation elements
  auto content_mat = (double *)(&fNavInd[fCurrent]);
  for (auto i = 0; i < 3; ++i)
    content_mat[i] = mat.Translation(i);
  for (auto i = 0; i < 9; ++i)
    content_mat[i + 3] = mat.Rotation(i);

  // Set new value for fCurrent
  fCurrent += 24;
  return new_mother;
}

NavIndex_t NavIndexTable::ValidateState(NavStatePath *state)
{
  // Decode the NavIndex_t
  unsigned char level            = state->GetLevel();
  int dind                       = 0;
  NavIndex_t nav_ind             = fWorld;
  VPlacedVolume const *pdaughter = nullptr;
  for (int i = 1; i < level + 1; ++i) {
    pdaughter = state->At(i);
    dind      = pdaughter->GetChildId();
    if (dind < 0) {
      std::runtime_error("=== EEE === Validate: incompatible daughter pointer");
      state->Print();
      return 0;
    }
    nav_ind = NavStateIndex::PushImpl(nav_ind, pdaughter);
    // nav_ind = Push(nav_ind, dind);
    assert(nav_ind > 0);
  }

  // Check if the physical volume is correct
  if (NavStateIndex::TopImpl(nav_ind) != state->Top()) {
    std::runtime_error("=== EEE === Validate: Top placed volume pointer mismatch");
    state->Print();
    return 0;
  }

  // Check if the current level is valid
  if (level != NavStateIndex::GetLevelImpl(nav_ind)) {
    std::runtime_error("=== EEE === Validate: Level mismatch");
    state->Print();
    return 0;
  }

  // Check if mother navigation index is consistent
  if (level > 0 && nav_ind != NavStateIndex::PushImpl(NavStateIndex::PopImpl(nav_ind), pdaughter)) {
    std::runtime_error("=== EEE === Validate: Navigation index inconsistency for Push/Pop");
    state->Print();
    return 0;
  }

  // Check if the number of daughters is correct
  if (NavStateIndex::GetNdaughtersImpl(nav_ind) != state->Top()->GetDaughters().size()) {
    std::runtime_error("=== EEE === Validate: Number of daughters mismatch");
    state->Print();
    return 0;
  }

  Transformation3D trans, trans_nav_ind;
  state->TopMatrix(trans);
  NavStateIndex::TopMatrixImpl(nav_ind, trans_nav_ind);
  if (!trans.operator==(trans_nav_ind)) {
    std::runtime_error("=== EEE === Validate: Transformation matrix mismatch");
    state->Print();
    return 0;
  }

  // success
  return nav_ind;
}

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom
