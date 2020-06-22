# Script to generate translation and rotation specializations of placed volumes,

rotation = [0x1B1, 0x18E, 0x076, 0x16A, 0x155, 0x0AD, 0x0DC, 0x0E3, 0x11B,
            0x0A1, 0x10A, 0x046, 0x062, 0x054, 0x111, 0x200]
translation = ["translation::kGeneric", "translation::kIdentity"]

header_string = """\
/// @file TransformationSpecializations.icc
/// @author Script generated.

#ifndef VECGEOM_NO_SPECIALIZATION
"""

specialization_string = """\
  if (trans_code == {:s} && rot_code == {:#05x}) {{
    return VolumeType::template Create<{:s}, {:#05x}>(
      logical_volume, transformation,
#ifdef VECCORE_CUDA
           id, copy_no, child_id,
#endif
      placement);
  }}
"""

generic_string = """\
#endif // No specialization

  return VolumeType::template Create<translation::kGeneric, rotation::kGeneric>(
      logical_volume, transformation,
#ifdef VECCORE_CUDA
           id, copy_no, child_id,
#endif
      placement);\
"""

output_file = open("../VecGeom/management/TransformationSpecializations.icc", "w")

output_file.write(header_string)
for r in rotation:
  for t in translation:
    output_file.write(specialization_string.format(t, r, t, r))
output_file.write(generic_string)
output_file.close()
