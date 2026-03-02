#ifndef _OSH_PARTICLE_PDG
#define _OSH_PARTICLE_PDG

// clang-format off
#define OSH_PART_PDG_NEUTRON     2112    /* neutron */
#define OSH_PART_PDG_PROTON      2212    /* proton */
#define OSH_PART_PDG_PIMINUS     -211    /* pi- */
#define OSH_PART_PDG_PIPLUS      211     /* pi+ */
#define OSH_PART_PDG_PIZERO      111     /* pi0 */
#define OSH_PART_PDG_ANEUTRON    -2112    /* antineutron */
#define OSH_PART_PDG_APROTON     -2212    /* antiproton */

#define OSH_PART_PDG_KMINUS      -321    /* K- */
#define OSH_PART_PDG_KPLUS       321     /* K+ */

#define OSH_PART_PDG_KZERO       311     /* K0 */
#define OSH_PART_PDG_KTLD        130 /* K~ */      // TODO check what this is, using K0_L now
#define OSH_PART_PDG_GAMMA       22      /* Gamma ray */
#define OSH_PART_PDG_ELECTRON    11      /* Electron */
#define OSH_PART_PDG_POSITRON    -11     /* Position */
#define OSH_PART_PDG_MUMINUS     13      /* Muon mu- */
#define OSH_PART_PDG_MUPLUS      -13     /* Muon mu+ */
#define OSH_PART_PDG_ENU         12      /* e-neutrino  */
#define OSH_PART_PDG_EANU        -12     /* e-antineutrino  */
#define OSH_PART_PDG_MUNU        14      /* mu-neutrino */
#define OSH_PART_PDG_MUANU       -14     /* mu-antineutrino */

#define OSH_PART_PDG_DEUTERON    1000010020  /* deuteron */
#define OSH_PART_PDG_TRITON      1000010030  /* triton */
#define OSH_PART_PDG_HE3         1000020030  /* He-3 */
#define OSH_PART_PDG_HE4         1000020040  /* He-4 */
#define OSH_PART_PDG_HIBASE      1000000000  /* base number to build hadron in ground state */

/* PDG codes 81–100 are reserved for generator-specific pseudoparticles and concepts. */
#define OSH_PART_PDG_NONE        81       /* zero-physics virtual particle for testing, currently also used if not set */
#define OSH_PART_PDG_INVALID     100      /* for error handling */
// clang-format on

#endif /* !_OSH_PARTICLE_PDG */