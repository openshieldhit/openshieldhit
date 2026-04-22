#ifndef _OSH_GEMCA2_VOXEL_HU
#define _OSH_GEMCA2_VOXEL_HU

#include <stdint.h>

#include "openshieldhit/status.h"

struct osh_material_workspace;

float osh_gemca_voxel_hu2rho(int16_t hu, char alg);
int osh_gemca_voxel_hu2idx(int16_t hu);
float osh_gemca_voxel_hu2wepl(int16_t hu, char alg);

enum osh_status osh_gemca_voxel_register_schneider_materials(struct osh_material_workspace *wm);
void osh_gemca_voxel_build_hu_lut(uint8_t lut[2601]);
enum osh_status osh_gemca_voxel_register_permatassari_materials(struct osh_material_workspace *wm);
void osh_gemca_voxel_build_hu_lut_permatassari2020(uint8_t lut[2601]);
float osh_gemca_voxel_hu2rho_permatassari2020(int16_t hu, int bin);

#endif /* _OSH_GEMCA2_VOXEL_HU */
