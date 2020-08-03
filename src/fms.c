/*
 Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
 the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
 reserved. See files LICENSE and NOTICE for details.

 This file is part of CEED, a collection of benchmarks, miniapps, software
 libraries and APIs for efficient high-order finite element and spectral
 element discretizations for exascale applications. For more information and
 source code availability see http://github.com/ceed.

 The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
 a collaborative effort of two U.S. Department of Energy organizations (Office
 of Science and the National Nuclear Security Administration) responsible for
 the planning and preparation of a capable exascale ecosystem, including
 software, applications, hardware, advanced system engineering and early
 testbed platforms, in support of the nation's exascale computing imperative.
*/

// We use strdup() - this is one way to get it:
#define _XOPEN_SOURCE 600

#include "fms.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


struct FmsMesh_private {
  FmsInt num_domains;    ///< Number of entries in the #domains array.
  FmsInt num_comps;      ///< Number of entries in the #comps array.
  FmsInt num_tags;       ///< Number of entries in the #tags array.
  FmsInt partition_id;   ///< A partition id.
  FmsInt num_partitions; ///< Total number of partitions.
  FmsDomain    *domains;  ///< An array of FmsDomain%s.
  FmsComponent *comps;    ///< An array of FmsComponent%s.
  FmsTag       *tags;     ///< An array of FmsTag%s.

  FmsInt num_domain_names; /**< Number of entries in the #domain_names array,
                                i.e. number of unique domain names. */
  FmsInt *domain_offsets;  /**< Offsets into the #domains array to the beginning
                                of a sequence of domains that share the same
                                name; the size of this array is
                                #num_domain_names + 1. */
  char **domain_names;     ///< An array of domain names.
};

struct FmsDomain_private {
  /// Number of entities stored in the #entities and #orientstions arrays
  FmsInt num_entities[FMS_NUM_ENTITY_TYPES];
  /** @brief Capacity (allocated size) of the #entities and #orientstions
      arrays, in number of entities. */
  FmsInt entities_caps[FMS_NUM_ENTITY_TYPES];
  /// Types of the side ids used in the arrays #entities.
  FmsIntType side_ids_type[FMS_NUM_ENTITY_TYPES];
  /// Mesh entities described as tuples of ids of side entities.
  void *entities[FMS_NUM_ENTITY_TYPES];
  /// Orientations of the side entities of all entities.
  FmsOrientation *orientations[FMS_NUM_ENTITY_TYPES];
  /** Pointer to the domain name stored in FmsMesh_private::domain_names. Not
      owned. */
  char *name;
  /// Id of the domain.
  FmsInt id;
  /** Dimension of the domain, i.e. the largest dimension of an entity in the
      domain. */
  FmsInt dim;
};

// An internal structure used by FmsComponent_private.
struct _FmsPart_private {
  FmsDomain domain; // Associated domain. Not owned.
  FmsInt num_entities[FMS_NUM_ENTITY_TYPES];
  // Types of the entities ids used by the arrays in entities_ids.
  FmsIntType entities_ids_type[FMS_NUM_ENTITY_TYPES];
  void *entities_ids[FMS_NUM_ENTITY_TYPES];
  // Orientations of the main entities relative to the entity definitions in the
  // domain. Only the orientation for the main entities will be defined and
  // used.
  FmsOrientation *orientations[FMS_NUM_ENTITY_TYPES];
};

struct FmsComponent_private {
  char *name; ///< Component name. Owned.
  FmsInt id; ///< Component id: index into the components array of its mesh.
  FmsInt dim; ///< Dimension of the mesh component.
  FmsInt num_main_entities; ///< Total number of main entities in all parts.
  FmsInt num_parts;
  FmsInt num_relations;
  struct _FmsPart_private *parts;
  FmsField coordinates; ///< Not owned.
  FmsInt *relations; ///< An array of FmsComponent ids.
};

struct FmsTag_private {
  char *name; ///< Tag name. Owned.
  FmsComponent comp; ///< The associated mesh component. Not Owned.
  FmsInt num_tag_descriptions; ///< Number of tag descriptions.
  FmsIntType tag_type; ///< Type of the entires of the #tags array.
  void *tags; /**< Integer tags assciated with all main entries of the
                   associated mesh component, ordered part-by-part. */
  void *described_tags; ///< List of tags that have descriptions.
  char **tag_descriptions;
};

struct FmsMetaData_private {
  FmsMetaDataType md_type; ///< Type of the meta-data.
  union {
    FmsIntType int_type;
    FmsScalarType scalar_type;
  } sub_type; /// Subtype of the meta-data.
  char *name; ///< Name of the meta-data.
  FmsInt num_entries; ///< Number of entries in the data array.
  void *data; ///< Meta-data array.
};

// Internal type used by FmsFieldDescriptor
struct _FmsFieldDescriptor_FixedOrder {
  FmsFieldType field_type;
  FmsBasisType basis_type;
  FmsInt       order;
};

struct FmsFieldDescriptor_private {
  char *name;
  FmsComponent component; ///< Not Owned.
  FmsFieldDescriptorType descr_type;
  union {
    void *any;
    struct _FmsFieldDescriptor_FixedOrder *fixed_order;
  } descriptor;
  FmsInt num_dofs;
};

struct FmsField_private {
  char *name;
  FmsMetaData mdata;
  FmsFieldDescriptor fd; ///< Not owned.
  FmsInt num_vec_comp;
  FmsLayoutType layout;
  FmsScalarType scalar_type;
  void *data;
};

struct FmsDataCollection_private {
  char *name;
  FmsMetaData mdata;
  FmsMesh mesh;
  FmsInt num_fds;
  FmsFieldDescriptor *fds;
  FmsInt num_fields;
  FmsField *fields;
};

const size_t FmsIntTypeSize[FMS_NUM_INT_TYPES] = {
  1, 2, 4, 8, 1, 2, 4, 8
};

const char * const FmsIntTypeNames[FMS_NUM_INT_TYPES] = {
  "FMS_INT8",
  "FMS_INT16",
  "FMS_INT32",
  "FMS_INT64",
  "FMS_UINT8",
  "FMS_UINT16",
  "FMS_UINT32",
  "FMS_UINT64",
};

const char * const FmsEntityTypeNames[FMS_NUM_ENTITY_TYPES] = {
  "FMS_VERTEX",
  "FMS_EDGE",
  "FMS_TRIANGLE",
  "FMS_QUADRILATERAL",
  "FMS_TETRAHEDRON",
  "FMS_HEXAHEDRON",
  "FMS_WEDGE",
  "FMS_PYRAMID"
};

const FmsInt FmsEntityDim[FMS_NUM_ENTITY_TYPES] = {
  0, 1, 2, 2, 3, 3, 3, 3
};

const FmsInt FmsEntityNumSides[FMS_NUM_ENTITY_TYPES] = {
  0, 2, 3, 4, 4, 6, 5, 5
};

const FmsInt FmsEntityNumVerts[FMS_NUM_ENTITY_TYPES] = {
  1, 2, 3, 4, 4, 8, 6, 5
};

const size_t FmsScalarTypeSize[FMS_NUM_SCALAR_TYPES] = {
  sizeof(float), sizeof(double), 2*sizeof(float), 2*sizeof(double)
};

const char * const FmsScalarTypeNames[FMS_NUM_SCALAR_TYPES] = {
  "FMS_FLOAT",
  "FMS_DOUBLE",
  "FMS_COMPLEX_FLOAT",
  "FMS_COMPLEX_DOUBLE",
};

const char * const FmsMetaDataTypeNames[FMS_NUM_METADATA_TYPES] = {
  "FMS_INTEGER",
  "FMS_SCALAR",
  "FMS_STRING",
  "FMS_META_DATA",
};


/* -------------------------------------------------------------------------- */
/* Auxiliary macros and static functions */
/* -------------------------------------------------------------------------- */

#ifdef NDEBUG
#define E_RETURN(code) return (code)
#else
static void FmsErrorDebug(int err_code, const char *func, const char *file,
                          int line) {
  fprintf(stderr,
          "\n\nFMS error encountered!\n"
          "Error code: %d\n"
          "Function:   %s\n"
          "File:       %s:%d\n\n"
          "Aborting ...\n\n", err_code, func, file, line);
  abort();
}
#ifndef _MSC_VER
#define FMS_PRETTY_FUNCTION __PRETTY_FUNCTION__
#else // for Visual Studio C++
#define FMS_PRETTY_FUNCTION __FUNCSIG__
#endif
#define E_RETURN(code) \
  do { \
    if (code) { \
      FmsErrorDebug(code, FMS_PRETTY_FUNCTION, __FILE__, __LINE__); \
    } \
    return (code); \
  } while (0)
#endif

static int UpdateError(int old_err, int new_err) {
  return old_err ? old_err : new_err;
}

static inline FmsInt FmsIntMin(FmsInt a, FmsInt b) {
  return (a <= b) ? a : b;
}

static inline FmsInt FmsIntMax(FmsInt a, FmsInt b) {
  return (a >= b) ? a : b;
}

// Return m=2^k, k-integer, such that: m >= n > m/2
static inline FmsInt NextPow2(FmsInt n) {
  FmsInt m = 1;
  while (m < n) { m *= 2; }
  return m;
}

static void FmsAbortNotImplemented() {
  fprintf(stderr, "\n\nNot implemented! Aborting ...\n\n");
  abort();
}


static void FmsInternalError() {
  fprintf(stderr, "\n\nFMS internal error detected! Aborting ...\n\n");
  abort();
}

static inline int FmsCopyString(const char *str, char **str_copy_p) {
  if (str == NULL) { *str_copy_p = NULL; return 0; }
  char *str_copy = strdup(str);
  if (str_copy == NULL) { return 1; }
  *str_copy_p = str_copy;
  return 0;
}

// Macro used in FmsIntConvertCopy()
#define FMS_ICC_IL(src_t,src,dst_t,dst,size) \
  do { \
    for (FmsInt idx = 0; idx < (size); idx++) { \
      ((dst_t*)(dst))[idx] = (dst_t)(((src_t*)(src))[idx]); \
    } \
  } while (0)

// Macro used in FmsIntConvertCopy()
#define FMS_ICC_SW(src_type,src,dst_t,dst,size) \
  switch (src_type) { \
  case FMS_INT8:   FMS_ICC_IL(  int8_t,src,dst_t,dst,size); break; \
  case FMS_INT16:  FMS_ICC_IL( int16_t,src,dst_t,dst,size); break; \
  case FMS_INT32:  FMS_ICC_IL( int32_t,src,dst_t,dst,size); break; \
  case FMS_INT64:  FMS_ICC_IL( int64_t,src,dst_t,dst,size); break; \
  case FMS_UINT8:  FMS_ICC_IL( uint8_t,src,dst_t,dst,size); break; \
  case FMS_UINT16: FMS_ICC_IL(uint16_t,src,dst_t,dst,size); break; \
  case FMS_UINT32: FMS_ICC_IL(uint32_t,src,dst_t,dst,size); break; \
  case FMS_UINT64: FMS_ICC_IL(uint64_t,src,dst_t,dst,size); break; \
  default: FmsAbortNotImplemented(); \
  } \
  break

// Copy items from one integer array to another, converting the entries.
static int FmsIntConvertCopy(FmsIntType src_type, const void *src,
                             FmsIntType dst_type, void *dst, FmsInt size) {
  if (src_type >= FMS_NUM_INT_TYPES || dst_type >= FMS_NUM_INT_TYPES) {
    E_RETURN(1);
  }
  switch (dst_type) {
  case FMS_INT8:   FMS_ICC_SW(src_type,src,  int8_t,dst,size);
  case FMS_INT16:  FMS_ICC_SW(src_type,src, int16_t,dst,size);
  case FMS_INT32:  FMS_ICC_SW(src_type,src, int32_t,dst,size);
  case FMS_INT64:  FMS_ICC_SW(src_type,src, int64_t,dst,size);
  case FMS_UINT8:  FMS_ICC_SW(src_type,src, uint8_t,dst,size);
  case FMS_UINT16: FMS_ICC_SW(src_type,src,uint16_t,dst,size);
  case FMS_UINT32: FMS_ICC_SW(src_type,src,uint32_t,dst,size);
  case FMS_UINT64: FMS_ICC_SW(src_type,src,uint64_t,dst,size);
  default: FmsAbortNotImplemented();
  }
  return 0;
}

// Macro used in FmsIntPermuteConvertCopy()
#define FMS_IPCC_IL(src_t,src,dst_t,dst,tuple_size,perm,tot_size) \
  do { \
    for (FmsInt i = 0; i < (tot_size); i += (tuple_size)) { \
      src_t *loc_src = (src_t*)(src) + i; \
      dst_t *loc_dst = (dst_t*)(dst) + i; \
      for (FmsInt j = 0; j < (tuple_size); j++) { \
        loc_dst[j] = (dst_t)(loc_src[perm[j]]); \
      } \
    } \
  } while (0)

// Macro used in FmsIntPermuteConvertCopy()
#define FMS_IPCC_SW(src_type,src,dst_t,dst,ts,perm,size) \
  switch (src_type) { \
  case FMS_INT8:   FMS_IPCC_IL(  int8_t,src,dst_t,dst,ts,perm,size); break; \
  case FMS_INT16:  FMS_IPCC_IL( int16_t,src,dst_t,dst,ts,perm,size); break; \
  case FMS_INT32:  FMS_IPCC_IL( int32_t,src,dst_t,dst,ts,perm,size); break; \
  case FMS_INT64:  FMS_IPCC_IL( int64_t,src,dst_t,dst,ts,perm,size); break; \
  case FMS_UINT8:  FMS_IPCC_IL( uint8_t,src,dst_t,dst,ts,perm,size); break; \
  case FMS_UINT16: FMS_IPCC_IL(uint16_t,src,dst_t,dst,ts,perm,size); break; \
  case FMS_UINT32: FMS_IPCC_IL(uint32_t,src,dst_t,dst,ts,perm,size); break; \
  case FMS_UINT64: FMS_IPCC_IL(uint64_t,src,dst_t,dst,ts,perm,size); break; \
  default: FmsAbortNotImplemented(); \
  } \
  break

// Copy items from one integer array to another, converting the entries and
// permuting the entries within tuples of given size.
static int FmsIntPermuteConvertCopy(FmsIntType src_type, const void *src,
                                    FmsIntType dst_type, void *dst,
                                    const int *perm, FmsInt tuple_size,
                                    FmsInt tot_size) {
  if (src_type >= FMS_NUM_INT_TYPES || dst_type >= FMS_NUM_INT_TYPES) {
    E_RETURN(1);
  }
  if (tuple_size == 0 || tot_size % tuple_size != 0) { E_RETURN(2); }
  const FmsInt ts = tuple_size;
  switch (dst_type) {
  case FMS_INT8:   FMS_IPCC_SW(src_type,src,  int8_t,dst,ts,perm,tot_size);
  case FMS_INT16:  FMS_IPCC_SW(src_type,src, int16_t,dst,ts,perm,tot_size);
  case FMS_INT32:  FMS_IPCC_SW(src_type,src, int32_t,dst,ts,perm,tot_size);
  case FMS_INT64:  FMS_IPCC_SW(src_type,src, int64_t,dst,ts,perm,tot_size);
  case FMS_UINT8:  FMS_IPCC_SW(src_type,src, uint8_t,dst,ts,perm,tot_size);
  case FMS_UINT16: FMS_IPCC_SW(src_type,src,uint16_t,dst,ts,perm,tot_size);
  case FMS_UINT32: FMS_IPCC_SW(src_type,src,uint32_t,dst,ts,perm,tot_size);
  case FMS_UINT64: FMS_IPCC_SW(src_type,src,uint64_t,dst,ts,perm,tot_size);
  default: FmsAbortNotImplemented();
  }
  return 0;
}

// Macro used in FmsIntMapCopy()
#define FMS_IMC_IL(src_t,src,dst_t,dst,map,size) \
  do { \
    for (FmsInt idx = 0; idx < (size); idx++) { \
      ((dst_t*)(dst))[idx] = ((dst_t*)(map))[((src_t*)(src))[idx]]; \
    } \
  } while (0)

// Macro used in FmsIntMapCopy()
#define FMS_IMC_SW(src_type,src,dst_t,dst,map,size) \
  switch (src_type) { \
  case FMS_INT8:   FMS_IMC_IL(  int8_t,src,dst_t,dst,map,size); break; \
  case FMS_INT16:  FMS_IMC_IL( int16_t,src,dst_t,dst,map,size); break; \
  case FMS_INT32:  FMS_IMC_IL( int32_t,src,dst_t,dst,map,size); break; \
  case FMS_INT64:  FMS_IMC_IL( int64_t,src,dst_t,dst,map,size); break; \
  case FMS_UINT8:  FMS_IMC_IL( uint8_t,src,dst_t,dst,map,size); break; \
  case FMS_UINT16: FMS_IMC_IL(uint16_t,src,dst_t,dst,map,size); break; \
  case FMS_UINT32: FMS_IMC_IL(uint32_t,src,dst_t,dst,map,size); break; \
  case FMS_UINT64: FMS_IMC_IL(uint64_t,src,dst_t,dst,map,size); break; \
  default: FmsAbortNotImplemented(); \
  } \
  break

static int FmsIntMapCopy(FmsIntType src_type, const void *src,
                         FmsIntType dst_type, void *dst, const void *map,
                         FmsInt size) {
  if (src_type >= FMS_NUM_INT_TYPES || dst_type >= FMS_NUM_INT_TYPES) {
    E_RETURN(1);
  }
  switch (dst_type) {
  case FMS_INT8:   FMS_IMC_SW(src_type,src,  int8_t,dst,map,size);
  case FMS_INT16:  FMS_IMC_SW(src_type,src, int16_t,dst,map,size);
  case FMS_INT32:  FMS_IMC_SW(src_type,src, int32_t,dst,map,size);
  case FMS_INT64:  FMS_IMC_SW(src_type,src, int64_t,dst,map,size);
  case FMS_UINT8:  FMS_IMC_SW(src_type,src, uint8_t,dst,map,size);
  case FMS_UINT16: FMS_IMC_SW(src_type,src,uint16_t,dst,map,size);
  case FMS_UINT32: FMS_IMC_SW(src_type,src,uint32_t,dst,map,size);
  case FMS_UINT64: FMS_IMC_SW(src_type,src,uint64_t,dst,map,size);
  default: FmsAbortNotImplemented();
  }
  return 0;
}

// Macro used in FmsEntitiesToSides()
#define FMS_E2S_TEMPL(idx_t) \
  do { \
    const idx_t *ents_ = (idx_t*)ents + loc_side_begin; \
    const idx_t *sides_ = sides; \
    for (FmsInt ent_id = 0; ent_id < num_ents; ent_id++) { \
      const idx_t *ent_off = &ents_[ent_id*ent_sides]; \
      FmsInt *res_ent = &res_[ent_id*res_stride]; \
      for (FmsInt loc_off = 0; loc_off < num_loc_sides; loc_off++) { \
        const idx_t side_id = ent_off[loc_off]; \
        const idx_t *side = &sides_[side_id*side_sides]; \
        FmsInt *res_side = &res_ent[loc_off*side_sides]; \
        for (FmsInt ss_id = 0; ss_id < side_sides; ss_id++) { \
          res_side[ss_id] = side[ss_id]; \
        } \
      } \
    } \
  } while (0); break

// Helper function for converting entities given as tuples of side ids to
// tuples of tuples where the inner tuples are the tuples defining the sides of
// the starting entities.
static int FmsEntitiesToSides(FmsIntType idx_type, FmsInt ent_sides,
                              FmsInt loc_side_begin, FmsInt num_loc_sides,
                              const void *ents, FmsInt side_sides,
                              const void *sides, FmsInt res_stride,
                              FmsInt res_off, FmsInt *res, FmsInt num_ents) {
  if (idx_type >= FMS_NUM_INT_TYPES) { E_RETURN(1); }
  FmsInt *res_ = res + res_off;
  switch (idx_type) {
  case FMS_INT8:   FMS_E2S_TEMPL(int8_t);
  case FMS_INT16:  FMS_E2S_TEMPL(int16_t);
  case FMS_INT32:  FMS_E2S_TEMPL(int32_t);
  case FMS_INT64:  FMS_E2S_TEMPL(int64_t);
  case FMS_UINT8:  FMS_E2S_TEMPL(uint8_t);
  case FMS_UINT16: FMS_E2S_TEMPL(uint16_t);
  case FMS_UINT32: FMS_E2S_TEMPL(uint32_t);
  case FMS_UINT64: FMS_E2S_TEMPL(uint64_t);
  default: FmsAbortNotImplemented();
  }
  return 0;
}

static inline int FmsIntersectTwoEdges(FmsInt *edge1, FmsInt *edge2,
                                       FmsInt *v_common) {
  int num_common = 0;
  if (edge1[0] == edge2[0] || edge1[0] == edge2[1]) {
    *v_common = edge1[0];
    num_common++;
  }
  if (edge1[0] != edge1[1]) {
    if (edge1[1] == edge2[0] || edge1[1] == edge2[1]) {
      *v_common = edge1[1];
      num_common++;
    }
  }
  return num_common;
}

static int FmsOrientFaceSides(FmsInt num_sides, FmsInt num_faces,
                              FmsInt *face_verts,
                              FmsOrientation *orientations) {
  // The size of 'face_verts' must be 2 x 'num_sides' x 'num_faces'.
  // The size of 'orientations' must be 'num_sides' x 'num_faces'.

  // Find edge orientations that if applied to 'face_verts' lead to an array of
  // the form: (a,b,b,c,c,a), for triangles, or (a,b,b,c,c,d,d,a), for quads.
  // If one of the edge orientations cannot be uniquely determined, that edge in
  // the 'orientations' array will remain unmodified, e.g. as initialized by the
  // user.

  const int n = num_sides;
  for (FmsInt k = 0; k < num_faces; k++) {
    FmsInt nc[4], vi[4]; // at most 4 vertices/sides
    for (int i = 0; i < n; i++) {
      nc[i] = FmsIntersectTwoEdges(face_verts+2*((i+n-1)%n),
                                   face_verts+2*i, vi+i);
      if (nc[i] == 0) { E_RETURN(1); }
    }
    for (int ei = 0; ei < n; ei++) {
      int ej = (ei+1)%n;
      if (face_verts[2*ei+0] == face_verts[2*ei+1] ||
          (nc[ei] != 1 && nc[ej] != 1)) {
        // Keep the value of 'orientations[ei]' as is, e.g. as defined by the
        // user.
      } else {
        if (nc[ei] == 1) {
          orientations[ei] = (vi[ei] == face_verts[2*ei+0]) ? 0 : 1;
        } else {
          orientations[ei] = (vi[ej] == face_verts[2*ei+1]) ? 0 : 1;
        }
      }
    }
    face_verts += 2*n;
    orientations += n;
  }
  return 0;
}

static inline void Sort3(const FmsInt *a, FmsInt *s) {
  FmsInt a0 = a[0], a1 = a[1], a2 = a[2];
  if (a0 > a1) {
    if (a2 <= a1) {
      s[0] = a2;
      s[1] = a1;
      s[2] = a0;
    } else {
      s[0] = a1;
      if (a2 < a0) {
        s[1] = a2;
        s[2] = a0;
      } else {
        s[1] = a0;
        s[2] = a2;
      }
    }
  } else { // a0 <= a1
    if (a1 <= a2) { // already sorted
      s[0] = a0;
      s[1] = a1;
      s[2] = a2;
    } else {
      if (a0 <= a2) {
        s[0] = a0;
        s[1] = a2;
      } else {
        s[0] = a2;
        s[1] = a0;
      }
      s[2] = a1;
    }
  }
}

static inline void Sort4(const FmsInt *a, FmsInt *s) {
  FmsInt b0, b1, b2, b3;
  {
    FmsInt a0 = a[0], a1 = a[1];
    if (a0 <= a1) {
      b0 = a0;
      b1 = a1;
    } else {
      b0 = a1;
      b1 = a0;
    }
    a0 = a[2];
    a1 = a[3];
    if (a0 <= a1) {
      b2 = a0;
      b3 = a1;
    } else {
      b2 = a1;
      b3 = a0;
    }
  }
  if (b0 <= b2) { // b0 <= b2 <= b3, and b0 <= b1
    s[0] = b0;
    if (b1 <= b2) {
      s[1] = b1;
      s[2] = b2;
      s[3] = b3;
    } else {
      s[1] = b2;
      if (b1 <= b3) { // b0 <= b2 < b1 <= b3
        s[2] = b1;
        s[3] = b3;
      } else { // b0 <= b2 <= b3 < b1
        s[2] = b3;
        s[3] = b1;
      }
    }
  } else { // b2 < b0 <= b1, and b2 <= b3
    s[0] = b2;
    if (b3 <= b0) { // b2 <= b3 <= b0 <= b1
      s[1] = b3;
      s[2] = b0;
      s[3] = b1;
    } else {
      s[1] = b0;
      if (b3 <= b1) { // b2 < b0 < b3 <= b1
        s[2] = b3;
        s[3] = b1;
      } else { // b2 < b0 <= b1 < b3
        s[2] = b1;
        s[3] = b3;
      }
    }
  }
}

// Assuming 'face1' and 'face2' are sorted
static inline int FmsIntersectSorted(int n1, int n2, FmsInt *face1,
                                     FmsInt *face2, FmsInt *edge_common) {
  int num_common = 0;
  FmsInt e1 = face1[0], e2 = face2[0];
  for (int i1 = 0, i2 = 0; 1; ) {
    if (e1 == e2) {
      num_common++;
      *edge_common = e1;
      do {
        if (++i2 == n2) { return num_common; }
        e2 = face2[i2];
      } while (e1 == e2);
      do {
        if (++i1 == n1) { return num_common; }
        e1 = face1[i1];
      } while (e1 < e2);
    } else if (e1 < e2) {
      do {
        if (++i1 == n1) { return num_common; }
        e1 = face1[i1];
      } while (e1 < e2);
    } else { // e1 > e2
      do {
        if (++i2 == n2) { return num_common; }
        e2 = face2[i2];
      } while (e1 > e2);
    }
  }
  // not reachable
}

static inline FmsOrientation GetTriangleOrientation(int *ne, FmsInt *ee,
    FmsInt *be) {
  // For orientation definitions, see the Doxygen documentation at the
  // definition of type FmsOrientation in fms.h.

  // To determine the orientation of a triangle we need at least two of the
  // edges to be uniquely determined, i.e. their ne[i] to be 1.

  // Orientations will be computed relative to (be[0],be[1],be[2]).

  // If there are repetitions in 'be', the orientation is unknown.
  if (be[0] == be[1] || be[0] == be[2] || be[1] == be[2]) {
    return FMS_ORIENTATION_UNKNOWN;
  }
  if (ne[0] == 1) {
    int i0 = 0;
    while (ee[0] != be[i0]) {
      i0++;
#ifndef NDEBUG
      if (i0 == 3) { FmsInternalError(); }
#endif
    }
    if (ne[1] == 1) {
      // Find the orientation of (ee[0],ee[1],*)
      return (be[(i0+1)%3] == ee[1]) ? 2*i0 : 2*i0+1;
    } else if (ne[2] == 1) {
      // Find the orientation of (ee[0],*,ee[2])
      return (be[(i0+2)%3] == ee[2]) ? 2*i0 : 2*i0+1;
    }
  } else if (ne[1] == 1 && ne[2] == 1) {
    // Find the orientation of (*,ee[1],ee[2])
    int i1 = 0;
    while (ee[1] != be[i1]) {
      i1++;
#ifndef NDEBUG
      if (i1 == 3) { FmsInternalError(); }
#endif
    }
    if (be[(i1+1)%3] == ee[2]) {
      static const FmsOrientation oo[3] = { 4, 0, 2 };
      return oo[i1];
    } else {
      static const FmsOrientation oo[3] = { 3, 5, 1 };
      return oo[i1];
    }
  }
  return FMS_ORIENTATION_UNKNOWN;
}

static inline FmsOrientation GetQuadOrientation(int *ne, FmsInt *ee,
    FmsInt *be, FmsInt *se) {
  // For orientation definitions, see the Doxygen documentation at the
  // definition of type FmsOrientation in fms.h.

  // To determine the orientation of a quad we need at least two consecutive of
  // its edges to be uniquely determined, i.e. their ne[i] to be 1.

  // Orientations will be computed relative to (be[0],be[1],be[2],be[3]).

  // If there are repeated edge indices in 'be' (or in its sorted version 'se')
  // then treat the orientation as unknown, even though it still may be possible
  // to determine the orientation.
  if (se[0] == se[1] || se[1] == se[2] || se[2] == se[3]) {
    return FMS_ORIENTATION_UNKNOWN;
  }
  // For now, if not all edges are uniquely deterined, i.e. ne[i] is not 1 for
  // some i, then treat the orientation as unknown, even though it still may be
  // possible to determine the orientation.
  if (ne[0] != 1 || ne[1] != 1 || ne[2] != 1 || ne[3] != 1) {
    return FMS_ORIENTATION_UNKNOWN;
  }
  int i0 = 0;
  while (ee[0] != be[i0]) {
    i0++;
#ifndef NDEBUG
    if (i0 == 4) { FmsInternalError(); }
#endif
  }
  // Find the orientation of (ee[0],ee[1],*,*)
  return (be[(i0+1)%4] == ee[1]) ? 2*i0 : 2*i0+1;
}

static int FmsOrientTetSides(FmsInt num_tets, FmsInt *tet_edges,
                             FmsOrientation *orientations) {
  for (FmsInt i = 0; i < num_tets; i++) {
    FmsInt tet_edges_sorted[4][3];
    for (int loc_face_id = 0; loc_face_id < 4; loc_face_id++) {
      Sort3(tet_edges+3*loc_face_id, tet_edges_sorted[loc_face_id]);
    }
    int nc[6];
    FmsInt ei[6];

    // Array describing the tet edges as intersections of two tet faces.
    // For the definitions of the faces and edges, see the diagram in fms.h.
    static const int tet_edge_face[6][2] = {
      // a=A^B   b=A^D   c=A^C   d=B^C   e=B^D   f=C^D
      // A=0, B=1, C=2, D=3
      {0, 1}, {0, 3}, {0, 2}, {1, 2}, {1, 3}, {2, 3}
    };
    // Array describing the tet faces through their edges.
    // For the definitions of the faces and edges, see the diagram in fms.h.
    static const int tet_face_edge[4][3] = {
      // a=0, b=1, c=2, d=3, e=4, f=5
      {0, 2, 1}, {0, 4, 3}, {2, 3, 5}, {1, 5, 4}
    };

    // Intersect the faces to find the edges
    for (int i = 0; i < 6; i++) {
      const int *f = tet_edge_face[i];
      nc[i] = FmsIntersectSorted(3, 3, tet_edges_sorted[f[0]],
                                 tet_edges_sorted[f[1]], ei+i);
      if (nc[i] == 0) { E_RETURN(1); }
    }

    // Determine side orientations
    for (int loc_face_id = 0; loc_face_id < 4; loc_face_id++) {
      const int *e = tet_face_edge[loc_face_id];
      int ne[3] = { nc[e[0]], nc[e[1]], nc[e[2]] };
      FmsInt ee[3] = { ei[e[0]], ei[e[1]], ei[e[2]] };
      FmsInt *be = tet_edges+3*loc_face_id;
      // Orientations will be computed relative to (be[0],be[1],be[2])
      FmsOrientation ori = GetTriangleOrientation(ne, ee, be);
      if (ori != FMS_ORIENTATION_UNKNOWN) {
        orientations[loc_face_id] = ori;
      } else {
        // The orientation of this face is unknown, so leave the value of
        // 'orientations[loc_face_id]' unmodified, e.g. as set by the user.
      }
    }
    tet_edges += 12;
    orientations += 4;
  }
  return 0;
}

static int FmsOrientHexSides(FmsInt num_hexes, FmsInt *hex_edges,
                             FmsOrientation *orientations) {
  for (FmsInt i = 0; i < num_hexes; i++) {
    FmsInt hex_edges_sorted[6][4];
    for (int loc_face_id = 0; loc_face_id < 6; loc_face_id++) {
      Sort4(hex_edges+4*loc_face_id, hex_edges_sorted[loc_face_id]);
    }
    int nc[12];
    FmsInt ei[12];

    // Array describing the hex edges as intersections of two hex faces.
    // For the definitions of the faces and edges, see the diagram in fms.h.
    static const int hex_edge_face[12][2] = {
      // a=A^C   b=A^F   c=A^D   d=A^E   e=B^C   f=B^F
      // g=B^D   h=B^E   i=C^E   j=C^F   k=D^F   l=D^E
      // A=0, B=1, C=2, D=3, E=4, F=5
      {0, 2}, {0, 5}, {0, 3}, {0, 4}, {1, 2}, {1, 5},
      {1, 3}, {1, 4}, {2, 4}, {2, 5}, {3, 5}, {3, 4}
    };
    // Array describing the hex faces through their edges.
    // For the definitions of the faces and edges, see the diagram in fms.h.
    static const int hex_face_edge[6][4] = {
      // a=0, b=1, c=2, d=3, e=4, f=5, g=6, h=7, i=8, j=9, k=10, l=11
      {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 9, 4, 8},
      {2, 11, 6, 10}, {3, 8, 7, 11}, {1, 10, 5, 9}
    };

    // Intersect the faces to find the edges
    for (int i = 0; i < 12; i++) {
      const int *f = hex_edge_face[i];
      nc[i] = FmsIntersectSorted(4, 4, hex_edges_sorted[f[0]],
                                 hex_edges_sorted[f[1]], ei+i);
      if (nc[i] == 0) { E_RETURN(1); }
    }

    // Determine side orientations
    for (int loc_face_id = 0; loc_face_id < 6; loc_face_id++) {
      const int *e = hex_face_edge[loc_face_id];
      int ne[4] = { nc[e[0]], nc[e[1]], nc[e[2]], nc[e[3]] };
      FmsInt ee[4] = { ei[e[0]], ei[e[1]], ei[e[2]], ei[e[3]] };
      FmsInt *be = hex_edges+4*loc_face_id;
      // Orientations will be computed relative to (be[0],be[1],be[2],be[3])
      FmsOrientation ori =
        GetQuadOrientation(ne, ee, be, hex_edges_sorted[loc_face_id]);
      if (ori != FMS_ORIENTATION_UNKNOWN) {
        orientations[loc_face_id] = ori;
      } else {
        // The orientation of this face is unknown, so leave the value of
        // 'orientations[loc_face_id]' unmodified, e.g. as set by the user.
      }
    }
    hex_edges += 24;
    orientations += 6;
  }
  return 0;
}

// Permute the edges/vertices of a triangle according to a given orientation.
// Note that the permutation arrays for the edges and the vertices are the same.
static inline void FmsPermuteTriBdr(FmsOrientation ori, const FmsInt *e_in,
                                    FmsInt *e_out) {
  static const int tri_perm[6][3] = {
    {0, 1, 2}, {0, 2, 1}, {1, 2, 0}, {1, 0, 2}, {2, 0, 1}, {2, 1, 0}
  };
  const int *perm = tri_perm[ori];
  for (int i = 0; i < 3; i++) {
    e_out[i] = e_in[perm[i]];
  }
}

// Permute the edges/vertices of a quad according to a given orientation.
// Note that the permutation arrays for the edges and the vertices are the same.
static inline void FmsPermuteQuadBdr(FmsOrientation ori, const FmsInt *e_in,
                                     FmsInt *e_out) {
  static const int quad_perm[8][4] = {
    {0, 1, 2, 3}, {0, 3, 2, 1}, {1, 2, 3, 0}, {1, 0, 3, 2},
    {2, 3, 0, 1}, {2, 1, 0, 3}, {3, 0, 1, 2}, {3, 2, 1, 0}
  };
  const int *perm = quad_perm[ori];
  for (int i = 0; i < 4; i++) {
    e_out[i] = e_in[perm[i]];
  }
}


/* -------------------------------------------------------------------------- */
/* Fms general functions */
/* -------------------------------------------------------------------------- */

int FmsGetInterfaceVersion(FmsInt *version) {
  if (!version) { E_RETURN(1); }
  *version = FMS_INTERFACE_VERSION;
  return 0;
}

int FmsGetIntTypeFromName(const char * const name, FmsIntType *type) {
  if(!name) E_RETURN(1);
  if(!type) E_RETURN(2);

  if(strcmp(FmsIntTypeNames[FMS_INT8], name) == 0)
    *type = FMS_INT8;
  else if(strcmp(FmsIntTypeNames[FMS_INT16], name) == 0)
    *type = FMS_INT16;
  else if(strcmp(FmsIntTypeNames[FMS_INT32], name) == 0)
    *type = FMS_INT32;
  else if(strcmp(FmsIntTypeNames[FMS_INT64], name) == 0)
    *type = FMS_INT64;
  else if(strcmp(FmsIntTypeNames[FMS_UINT8], name) == 0)
    *type = FMS_UINT8;
  else if(strcmp(FmsIntTypeNames[FMS_UINT16], name) == 0)
    *type = FMS_UINT16;
  else if(strcmp(FmsIntTypeNames[FMS_UINT32], name) == 0)
    *type = FMS_UINT32;
  else if(strcmp(FmsIntTypeNames[FMS_UINT64], name) == 0)
    *type = FMS_UINT64;
  else
    E_RETURN(3);

  return 0;
}

int FmsGetScalarTypeFromName(const char * const name, FmsScalarType *type) {
  if(!name) E_RETURN(1);
  if(!type) E_RETURN(2);

  if(strcmp(FmsScalarTypeNames[FMS_FLOAT], name) == 0)
    *type = FMS_FLOAT;
  else if(strcmp(FmsScalarTypeNames[FMS_DOUBLE], name) == 0)
    *type = FMS_DOUBLE;
  else if(strcmp(FmsScalarTypeNames[FMS_COMPLEX_FLOAT], name) == 0)
    *type = FMS_COMPLEX_FLOAT;
  else if(strcmp(FmsScalarTypeNames[FMS_COMPLEX_DOUBLE], name) == 0)
    *type = FMS_COMPLEX_DOUBLE;
  else
    E_RETURN(3);

  return 0;
}

int FmsGetEntityTypeFromName(const char * const name, FmsEntityType *ent_type) {
  if(!name) E_RETURN(1);
  if(!ent_type) E_RETURN(2);

  if(strcmp(FmsEntityTypeNames[FMS_VERTEX], name) == 0)
    *ent_type = FMS_VERTEX;
  else if(strcmp(FmsEntityTypeNames[FMS_EDGE], name) == 0)
    *ent_type = FMS_EDGE;
  else if(strcmp(FmsEntityTypeNames[FMS_TRIANGLE], name) == 0)
    *ent_type = FMS_TRIANGLE;
  else if(strcmp(FmsEntityTypeNames[FMS_QUADRILATERAL], name) == 0)
    *ent_type = FMS_QUADRILATERAL;
  else if(strcmp(FmsEntityTypeNames[FMS_TETRAHEDRON], name) == 0)
    *ent_type = FMS_TETRAHEDRON;
  else if(strcmp(FmsEntityTypeNames[FMS_HEXAHEDRON], name) == 0)
    *ent_type = FMS_HEXAHEDRON;
  else if(strcmp(FmsEntityTypeNames[FMS_WEDGE], name) == 0)
    *ent_type = FMS_WEDGE;
  else if(strcmp(FmsEntityTypeNames[FMS_PYRAMID], name) == 0)
    *ent_type = FMS_PYRAMID;
  else
    E_RETURN(3);

  return 0;
}

int FmsGetMetaDataTypeFromName(const char * const name, FmsMetaDataType *type) {
  if(!name) E_RETURN(1);
  if(!type) E_RETURN(2);

  if(strcmp(FmsMetaDataTypeNames[FMS_INTEGER], name) == 0)
    *type = FMS_INTEGER;
  else if(strcmp(FmsMetaDataTypeNames[FMS_SCALAR], name) == 0)
    *type = FMS_SCALAR;
  else if(strcmp(FmsMetaDataTypeNames[FMS_STRING], name) == 0)
    *type = FMS_STRING;
  else if(strcmp(FmsMetaDataTypeNames[FMS_META_DATA], name) == 0)
    *type = FMS_META_DATA;
  else
    E_RETURN(3);
    
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsMesh functions: construction */
/* -------------------------------------------------------------------------- */

int FmsMeshConstruct(FmsMesh *mesh) {
  if (!mesh) { E_RETURN(1); }
  FmsMesh m;
  m = calloc(1, sizeof(*m));
  if (!m) { E_RETURN(2); }
  m->num_partitions = 1;
  *mesh = m;
  return 0;
}

int FmsMeshFinalize(FmsMesh mesh) {
  if (!mesh) { E_RETURN(1); }
  // Generate the side entity orientations for all entities in all domains
  for (FmsInt domain_id = 0; domain_id < mesh->num_domains; domain_id++) {
    FmsDomain domain = mesh->domains[domain_id];
    // Triangles
    if (domain->num_entities[FMS_TRIANGLE]) {
      const FmsEntityType et = FMS_TRIANGLE;
      const FmsIntType side_ids_type = domain->side_ids_type[et];
      if (side_ids_type != domain->side_ids_type[FMS_EDGE]) { E_RETURN(2); }
      const size_t sizeof_side_ids = FmsIntTypeSize[side_ids_type];
      const FmsInt num_tri = domain->num_entities[et];
      const char *all_tris = domain->entities[et];
      const void *all_edges = domain->entities[FMS_EDGE];
      FmsOrientation *oris = domain->orientations[et];
      const int bsize = 1024;
      FmsInt tri_verts[6*bsize];
      for (FmsInt off = 0; off < num_tri; off += bsize) {
        const FmsInt sz = FmsIntMin(bsize, num_tri-off);
        FmsEntitiesToSides(side_ids_type, 3, 0, 3, all_tris, 2, all_edges,
                           6, 0, tri_verts, sz);
        int err = FmsOrientFaceSides(3, sz, tri_verts, oris);
        if (err) { E_RETURN(err); }
        // TODO: create/update a list of the unknown side orientations
        all_tris += 3*sz*sizeof_side_ids;
        oris += 3*sz;
      }
    }
    // Quads
    if (domain->num_entities[FMS_QUADRILATERAL]) {
      const FmsEntityType et = FMS_QUADRILATERAL;
      const FmsIntType side_ids_type = domain->side_ids_type[et];
      if (side_ids_type != domain->side_ids_type[FMS_EDGE]) { E_RETURN(2); }
      const size_t sizeof_side_ids = FmsIntTypeSize[side_ids_type];
      const FmsInt num_quads = domain->num_entities[et];
      const char *all_quads = domain->entities[et];
      const void *all_edges = domain->entities[FMS_EDGE];
      FmsOrientation *oris = domain->orientations[et];
      const int bsize = 1024;
      FmsInt quad_verts[8*bsize];
      for (FmsInt off = 0; off < num_quads; off += bsize) {
        const FmsInt sz = FmsIntMin(bsize, num_quads-off);
        FmsEntitiesToSides(side_ids_type, 4, 0, 4, all_quads, 2, all_edges,
                           8, 0, quad_verts, sz);
        int err = FmsOrientFaceSides(4, sz, quad_verts, oris);
        if (err) { E_RETURN(err); }
        // TODO: create/update a list of the unknown side orientations
        all_quads += 4*sz*sizeof_side_ids;
        oris += 4*sz;
      }
    }
    // Tets
    if (domain->num_entities[FMS_TETRAHEDRON]) {
      const FmsEntityType et = FMS_TETRAHEDRON;
      const FmsIntType side_ids_type = domain->side_ids_type[et];
      if (side_ids_type != domain->side_ids_type[FMS_TRIANGLE]) { E_RETURN(2); }
      const size_t sizeof_side_ids = FmsIntTypeSize[side_ids_type];
      const FmsInt num_tets = domain->num_entities[et];
      const char *all_tets = domain->entities[et];
      const void *all_tris = domain->entities[FMS_TRIANGLE];
      FmsOrientation *oris = domain->orientations[et];
      const int bsize = 512;
      FmsInt tet_edges[12*bsize];
      for (FmsInt off = 0; off < num_tets; off += bsize) {
        const FmsInt sz = FmsIntMin(bsize, num_tets-off);
        FmsEntitiesToSides(side_ids_type, 4, 0, 4, all_tets, 3, all_tris,
                           12, 0, tet_edges, sz);
        int err = FmsOrientTetSides(sz, tet_edges, oris);
        if (err) { E_RETURN(err); }
        // TODO: create/update a list of the unknown side orientations
        all_tets += 4*sz*sizeof_side_ids;
        oris += 4*sz;
      }
    }
    // Hexes
    if (domain->num_entities[FMS_HEXAHEDRON]) {
      const FmsEntityType et = FMS_HEXAHEDRON;
      const FmsIntType side_ids_type = domain->side_ids_type[et];
      if (side_ids_type != domain->side_ids_type[FMS_QUADRILATERAL]) {
        E_RETURN(2);
      }
      const size_t sizeof_side_ids = FmsIntTypeSize[side_ids_type];
      const FmsInt num_hexes = domain->num_entities[et];
      const char *all_hexes = domain->entities[et];
      const void *all_quads = domain->entities[FMS_QUADRILATERAL];
      FmsOrientation *oris = domain->orientations[et];
      const int bsize = 256;
      FmsInt hex_edges[24*bsize];
      for (FmsInt off = 0; off < num_hexes; off += bsize) {
        const FmsInt sz = FmsIntMin(bsize, num_hexes-off);
        FmsEntitiesToSides(side_ids_type, 6, 0, 6, all_hexes, 4, all_quads,
                           24, 0, hex_edges, sz);
        int err = FmsOrientHexSides(sz, hex_edges, oris);
        if (err) { E_RETURN(err); }
        // TODO: create/update a list of the unknown side orientations
        all_hexes += 6*sz*sizeof_side_ids;
        oris += 6*sz;
      }
    }
    // Wedges
    if (domain->num_entities[FMS_WEDGE]) {
      // TODO ...
      FmsAbortNotImplemented();
    }
    // Pyramids
    if (domain->num_entities[FMS_PYRAMID]) {
      // TODO ...
      FmsAbortNotImplemented();
    }
  }
  // TODO ...
  return 0;
}

int FmsMeshValidate(FmsMesh mesh) {
  if (!mesh) { E_RETURN(1); }
  // Perform some validation
  if (mesh->domain_offsets[0] != 0) { E_RETURN(2); }
  if (mesh->domain_offsets[mesh->num_domain_names] != mesh->num_domains) {
    E_RETURN(3);
  }
  if (mesh->partition_id >= mesh->num_partitions) { E_RETURN(4); }
  // Check domains for consistency
  for (FmsInt name_id = 0; name_id < mesh->num_domain_names; name_id++) {
    const FmsInt name_start = mesh->domain_offsets[name_id];
    const FmsInt name_end = mesh->domain_offsets[name_id+1];
    for (FmsInt domain_id = name_start; domain_id < name_end; domain_id++) {
      FmsDomain domain = mesh->domains[domain_id];
      if (domain->name != mesh->domain_names[name_id]) { E_RETURN(5); }
      if (name_start + domain->id != domain_id) { E_RETURN(6); }
      // TODO: Check that all side ids are in the expected ranges
    }
  }
  // Check that all components use domains from the same mesh
  char **domain_names = mesh->domain_names;
  for (FmsInt comp_id = 0; comp_id < mesh->num_comps; comp_id++) {
    FmsComponent comp = mesh->comps[comp_id];
    if (comp->id != comp_id) { E_RETURN(7); }
    if (comp->dim == FMS_INVALID_DIM) { E_RETURN(8); }
    for (FmsInt part_id = 0; part_id < comp->num_parts; part_id++) {
      FmsDomain domain = comp->parts[part_id].domain;
      FmsInt domain_name_id = 0;
      for ( ; 1; domain_name_id++) {
        if (domain_name_id == mesh->num_domain_names) { E_RETURN(9); }
        if (strcmp(domain->name, domain_names[domain_name_id]) == 0) { break; }
      }
      FmsInt b_offset = mesh->domain_offsets[domain_name_id];
      FmsInt num_domains = mesh->domain_offsets[domain_name_id+1] - b_offset;
      if (domain->id >= num_domains) { E_RETURN(10); }
      if (mesh->domains[b_offset + domain->id] != domain) { E_RETURN(11); }
    }
  }
  // Check that all tags use components from the same mesh
  for (FmsInt tag_id = 0; tag_id < mesh->num_tags; tag_id++) {
    FmsTag tag = mesh->tags[tag_id];
    FmsInt comp_id = tag->comp->id;
    if (comp_id >= mesh->num_comps) { E_RETURN(12); }
    if (mesh->comps[comp_id] != tag->comp) { E_RETURN(13); }
  }
  return 0;
}

int FmsMeshDestroy(FmsMesh *mesh) {
  if (!mesh) { E_RETURN(1); }
  FmsMesh m = *mesh;
  if (m) {
    // Delete domains
    for (FmsInt i = 0; i < m->num_domains; i++) {
      FmsDomain domain = m->domains[i];
      for (FmsInt j = 0; j < FMS_NUM_ENTITY_TYPES; j++) {
        free(domain->orientations[j]);
        free(domain->entities[j]);
      }
      // domain->name is not owned
      free(domain);
    }
    free(m->domains);
    // Delete components
    for (FmsInt i = 0; i < m->num_comps; i++) {
      FmsComponent comp = m->comps[i];
      free(comp->relations);
      // comp->coordinates is not owned
      // Delete parts
      for (FmsInt j = 0; j < comp->num_parts; j++) {
        struct _FmsPart_private *part = &comp->parts[j];
        for (FmsInt k = 0; k < FMS_NUM_ENTITY_TYPES; k++) {
          free(part->orientations[k]);
          free(part->entities_ids[k]);
        }
      }
      free(comp->parts);
      free(comp->name);
      free(comp);
    }
    free(m->comps);
    // Delete tags
    for (FmsInt i = 0; i < m->num_tags; i++) {
      FmsTag tag = m->tags[i];
      for (FmsInt j = 0; j < tag->num_tag_descriptions; j++) {
        free(tag->tag_descriptions[j]);
      }
      free(tag->tag_descriptions);
      free(tag->described_tags);
      free(tag->tags);
      free(tag->name);
      free(tag);
    }
    free(m->tags);
    // Delete domain names and offsets
    for (FmsInt i = 0; i < m->num_domain_names; i++) {
      free(m->domain_names[i]);
    }
    free(m->domain_names);
    free(m->domain_offsets);
    // Delete the mesh struct
    free(m);
    *mesh = NULL;
  }
  return 0;
}

int FmsMeshSetPartitionId(FmsMesh mesh, FmsInt partition_id,
                          FmsInt num_partitions) {
  if (!mesh) { E_RETURN(1); }
  // if (partition_id < 0) { E_RETURN(2); } // if FmsInt is signed
  if (partition_id >= num_partitions) { E_RETURN(2); }
  mesh->partition_id = partition_id;
  mesh->num_partitions = num_partitions;
  return 0;
}

int FmsMeshAddDomains(FmsMesh mesh, const char *domain_name, FmsInt num_domains,
                      FmsDomain **domains) {
  if (!mesh) { E_RETURN(1); }
  if (!domain_name) { E_RETURN(2); }
  if (!domains) { E_RETURN(3); }
  FmsInt num_domain_names = mesh->num_domain_names;
  char **domain_names = mesh->domain_names;
  // Make sure the 'domain_name' is not already used
  for (FmsInt i = 0; i < num_domain_names; i++) {
    if (strcmp(domain_name, domain_names[i]) == 0) { E_RETURN(4); }
  }
  FmsInt *domain_offsets = mesh->domain_offsets;
  // Reallocate mesh->domain_names and mesh->domain_offsets, if necessary
  FmsInt nc = NextPow2(num_domain_names+1);
  if (num_domain_names <= nc/2) {
    domain_names = realloc(domain_names, nc*sizeof(*domain_names));
    if (domain_names == NULL) { E_RETURN(5); }
    mesh->domain_names = domain_names;
    domain_offsets = realloc(domain_offsets, (nc+1)*sizeof(FmsInt));
    if (domain_offsets == NULL) { E_RETURN(5); }
    domain_offsets[0] = 0;
    mesh->domain_offsets = domain_offsets;
  }
  FmsDomain *m_domains = mesh->domains;
  const FmsInt m_num_domains = mesh->num_domains;
  // Reallocate mesh->domains, if necessary
  if (num_domains != 0) {
    nc = NextPow2(m_num_domains + num_domains);
    if (m_num_domains <= nc/2) {
      m_domains = realloc(m_domains, nc*sizeof(FmsDomain));
      if (m_domains == NULL) { E_RETURN(5); }
      mesh->domains = m_domains;
    }
  }
  char *domain_name_copy;
  if (FmsCopyString(domain_name, &domain_name_copy)) { E_RETURN(5); }
  // On error, call free(domain_name_copy)
  for (FmsInt i = 0; i < num_domains; i++) {
    FmsDomain domain;
    domain = calloc(1, sizeof(*domain));
    if (domain == NULL) {
      while (i != 0) { free(m_domains[m_num_domains + (--i)]); }
      free(domain_name_copy);
      E_RETURN(5);
    }
    // Init the domain name, id, and dim
    domain->name = domain_name_copy;
    domain->id = i;
    domain->dim = FMS_INVALID_DIM;
    m_domains[m_num_domains + i] = domain;
  }
  domain_names[num_domain_names] = domain_name_copy;
  mesh->num_domain_names = ++num_domain_names;
  mesh->num_domains = domain_offsets[num_domain_names] = m_num_domains +
                      num_domains;
  *domains = m_domains + m_num_domains;
  return 0;
}

int FmsMeshAddComponent(FmsMesh mesh, const char *comp_name,
                        FmsComponent *comp) {
  if (!mesh) { E_RETURN(1); }
  if (!comp) { E_RETURN(2); }
  FmsInt num_comps = mesh->num_comps;
  FmsComponent *comps = mesh->comps;
  // Reallocate mesh->comps, if necessary
  FmsInt nc = NextPow2(num_comps+1);
  if (num_comps <= nc/2) {
    comps = realloc(comps, nc*sizeof(*comps));
    if (comps == NULL) { E_RETURN(3); }
    mesh->comps = comps;
  }
  FmsComponent component;
  component = calloc(1, sizeof(*component));
  if (component == NULL) { E_RETURN(3); }
  // On error, call free(component)
  if (FmsCopyString(comp_name, &component->name)) {
    free(component); E_RETURN(3);
  }
  component->id = num_comps;
  component->dim = FMS_INVALID_DIM;
  comps[num_comps] = component;
  mesh->num_comps = num_comps+1;
  *comp = component;
  return 0;
}

int FmsMeshAddTag(FmsMesh mesh, const char *tag_name, FmsTag *tag) {
  if (!mesh) { E_RETURN(1); }
  if (!tag) { E_RETURN(2); }
  FmsInt num_tags = mesh->num_tags;
  FmsTag *tags = mesh->tags;
  // Reallocate mesh->tags, if necessary
  FmsInt nc = NextPow2(num_tags+1);
  if (num_tags <= nc/2) {
    tags = realloc(tags, nc*sizeof(*tags));
    if (tags == NULL) { E_RETURN(3); }
    mesh->tags = tags;
  }
  FmsTag t;
  t = calloc(1, sizeof(*t));
  if (t == NULL) { E_RETURN(3); }
  // On error, call free(t)
  if (FmsCopyString(tag_name, &t->name)) { free(t); E_RETURN(3); }
  t->tag_type = FMS_NUM_INT_TYPES; // invalid type
  tags[num_tags] = t;
  mesh->num_tags = num_tags+1;
  *tag = t;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsDomain functions: construction */
/* -------------------------------------------------------------------------- */

int FmsDomainSetNumVertices(FmsDomain domain, FmsInt num_verts) {
  if (!domain) { E_RETURN(1); }
  domain->num_entities[FMS_VERTEX] = num_verts;
  if (domain->dim == FMS_INVALID_DIM) { domain->dim = 0; }
  return 0;
}

int FmsDomainSetNumEntities(FmsDomain domain, FmsEntityType type,
                            FmsIntType id_store_type, FmsInt num_ents) {
  if (!domain) { E_RETURN(1); }
  if (type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  if (type == FMS_VERTEX) {
    domain->num_entities[type] = num_ents;
    if (domain->dim == FMS_INVALID_DIM) { domain->dim = 0; }
    return 0;
  }
  if (id_store_type >= FMS_NUM_INT_TYPES) { E_RETURN(3); }
  const size_t sizeof_side_ids = FmsIntTypeSize[id_store_type];
  const FmsInt num_sides = FmsEntityNumSides[type];
  domain->num_entities[type] = 0;
  domain->entities_caps[type] = 0;
  free(domain->orientations[type]);
  free(domain->entities[type]);
  if (num_ents != 0) {
    FmsInt tot_sides = num_ents*num_sides;
    domain->entities[type] = malloc(tot_sides*sizeof_side_ids);
    if (domain->entities[type] == NULL) {
      domain->orientations[type] = NULL;
      E_RETURN(4);
    }
    domain->orientations[type] = malloc(tot_sides*sizeof(FmsOrientation));
    if (domain->orientations[type] == NULL) {
      free(domain->entities[type]);
      domain->entities[type] = NULL;
      E_RETURN(4);
    }
    FmsOrientation *orientations = domain->orientations[type];
    for (FmsInt i = 0; i < tot_sides; i++) {
      orientations[i] = FMS_ORIENTATION_UNKNOWN;
    }
  } else {
    domain->entities[type] = NULL;
    domain->orientations[type] = NULL;
  }
  domain->entities_caps[type] = num_ents;
  domain->side_ids_type[type] = id_store_type;
  if (domain->dim == FMS_INVALID_DIM) {
    domain->dim = FmsEntityDim[type];
  } else {
    domain->dim = FmsIntMax(domain->dim, FmsEntityDim[type]);
  }
  return 0;
}

int FmsDomainGetEntitiesArray(FmsDomain domain, FmsEntityType type,
                              void **ents) {
  if (!domain) { E_RETURN(1); }
  if (type == FMS_VERTEX || type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  if (ents == NULL) { E_RETURN(3); }
  // This function assumes that all domain->entities_caps[type] entities in the
  // returned array *ents will be set.
  domain->num_entities[type] = domain->entities_caps[type];
  *ents = domain->entities[type];
  return 0;
}

int FmsDomainGetOrientationsArray(FmsDomain domain, FmsEntityType type,
                                  FmsOrientation **side_orients) {
  if (!domain) { E_RETURN(1); }
  if (type == FMS_VERTEX || type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  if (side_orients == NULL) { E_RETURN(3); }
  *side_orients = domain->orientations[type];
  return 0;
}

int FmsDomainAddEntities(FmsDomain domain, FmsEntityType type,
                         FmsEntityReordering reordering, FmsIntType ent_id_type,
                         const void *ents, FmsInt num_ents) {
  if (!domain) { E_RETURN(1); }
  if (type == FMS_VERTEX || type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  const FmsInt num_entities = domain->num_entities[type];
  if (num_entities + num_ents > domain->entities_caps[type]) { E_RETURN(3); }
  const FmsIntType side_ids_type = domain->side_ids_type[type];
  const size_t sizeof_side_ids = FmsIntTypeSize[side_ids_type];
  const FmsInt num_sides = FmsEntityNumSides[type];
  void *ents_dst = (char*)(domain->entities[type]) +
                   num_entities*num_sides*sizeof_side_ids;
  const int *perm = reordering ? reordering[type] : NULL;
  if (perm == NULL) {
    FmsIntConvertCopy(ent_id_type, ents, side_ids_type, ents_dst,
                      num_ents*num_sides);
  } else {
    FmsIntPermuteConvertCopy(ent_id_type, ents, side_ids_type, ents_dst,
                             perm, num_sides, num_ents*num_sides);
  }
  domain->num_entities[type] = num_entities + num_ents;
  return 0;
}

int FmsDomainAddOrientation(FmsDomain domain, FmsEntityType type, FmsInt ent_id,
                            FmsInt loc_side_id, FmsOrientation side_orient) {
  if (!domain) { E_RETURN(1); }
  if (type == FMS_VERTEX || type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  if (ent_id >= domain->entities_caps[type]) { E_RETURN(3); }
  const FmsInt num_sides = FmsEntityNumSides[type];
  if (loc_side_id >= num_sides) { E_RETURN(4); }
  FmsOrientation *orientations = domain->orientations[type];
  orientations[ent_id*num_sides + loc_side_id] = side_orient;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsComponent functions: construction */
/* -------------------------------------------------------------------------- */

int FmsComponentAddDomain(FmsComponent comp, FmsDomain domain) {
  if (!comp) { E_RETURN(1); }
  if (!domain) { E_RETURN(2); }
  const FmsInt num_parts = comp->num_parts;
  const FmsInt comp_dim = comp->dim;
  const FmsInt dim = domain->dim;
  if (comp_dim != FMS_INVALID_DIM && comp_dim != dim) { E_RETURN(3); }
  struct _FmsPart_private *parts = comp->parts;
  // Reallocate comp->parts, if necessary
  FmsInt nc = NextPow2(num_parts + 1);
  if (num_parts <= nc/2) {
    parts = realloc(parts, nc*sizeof(*parts));
    if (parts == NULL) { E_RETURN(4); }
    comp->parts = parts;
  }
  struct _FmsPart_private *part = parts + num_parts;
  // Initialize part, counting the number of main entities
  memset(part, 0, sizeof(*part));
  part->domain = domain;
  FmsInt num_main_entities = comp->num_main_entities;
  for (FmsInt t = 0; t < FMS_NUM_ENTITY_TYPES; t++) {
    FmsInt num_entities = domain->num_entities[t];
    part->num_entities[t] = num_entities;
    if (FmsEntityDim[t] == dim) {
      num_main_entities += num_entities;
    }
    // the value of part->entities_ids_type[t] is not used
    // part->entities_ids[t] is already set to NULL
    // part->orientations[t] is already set to NULL
  }
  // Update comp
  if (comp_dim == FMS_INVALID_DIM) {
    comp->dim = dim;
  }
  comp->num_main_entities = num_main_entities;
  comp->num_parts = num_parts + 1;
  return 0;
}

int FmsComponentAddPart(FmsComponent comp, FmsDomain domain, FmsInt *part_id) {
  if (!comp) { E_RETURN(1); }
  if (!domain) { E_RETURN(2); }
  if (!part_id) { E_RETURN(3); }
  const FmsInt num_parts = comp->num_parts;
  struct _FmsPart_private *parts = comp->parts;
  // Reallocate comp->parts, if necessary
  FmsInt nc = NextPow2(num_parts + 1);
  if (num_parts <= nc/2) {
    parts = realloc(parts, nc*sizeof(*parts));
    if (parts == NULL) { E_RETURN(4); }
    comp->parts = parts;
  }
  struct _FmsPart_private *part = parts + num_parts;
  // Initialize part
  memset(part, 0, sizeof(*part));
  part->domain = domain;
  comp->num_parts = num_parts + 1;
  *part_id = num_parts;
  return 0;
}

int FmsComponentAddPartEntities(FmsComponent comp, FmsInt part_id,
                                FmsEntityType type, FmsIntType id_store_type,
                                FmsIntType ent_id_type, FmsIntType orient_type,
                                const FmsOrientation *inv_orient_map,
                                const void *ents, const void *ent_orients,
                                FmsInt num_ents) {
  if (!comp) { E_RETURN(1); }
  if (part_id >= comp->num_parts) { E_RETURN(2); }
  if (type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(3); }
  if (comp->dim != FMS_INVALID_DIM && comp->dim != FmsEntityDim[type]) {
    E_RETURN(3);
  }
  if (id_store_type >= FMS_NUM_INT_TYPES) { E_RETURN(4); }
  if (ent_id_type >= FMS_NUM_INT_TYPES) { E_RETURN(5); }
  if (orient_type >= FMS_NUM_INT_TYPES) { E_RETURN(6); }
  struct _FmsPart_private *part = comp->parts + part_id;
  const FmsInt num_entities = part->num_entities[type];
  void *entities_ids = part->entities_ids[type];
  FmsOrientation *orientations = part->orientations[type];
  if (num_entities != 0) {
    if (entities_ids == NULL) { E_RETURN(7); }
    if (orientations == NULL) { E_RETURN(8); }
    if (ents == NULL) { E_RETURN(9); }
    if (ent_orients == NULL) { E_RETURN(10); }
    if (part->entities_ids_type[type] != id_store_type) { E_RETURN(11); }
  } else if (ents == NULL) {
    if (num_ents != part->domain->num_entities[type]) { E_RETURN(12); }
  }
  const size_t sizeof_ent_ids = FmsIntTypeSize[id_store_type];
  if (num_ents != 0) {
    const FmsInt nc = NextPow2(num_entities + num_ents);
    if (ents != NULL) {
      if (num_entities <= nc/2) {
        entities_ids = realloc(entities_ids, nc*sizeof_ent_ids);
        if (entities_ids == NULL) { E_RETURN(13); }
        part->entities_ids[type] = entities_ids;
      }
      FmsIntConvertCopy(ent_id_type, ents, id_store_type,
                        (char*)entities_ids + num_entities*sizeof_ent_ids,
                        num_ents);
    } else {
      // part->entities_ids[type] is already NULL since num_entities is 0
    }
    if (ent_orients != NULL) {
      if (num_entities <= nc/2) {
        orientations = realloc(orientations, nc*sizeof(FmsOrientation));
        if (orientations == NULL) { E_RETURN(13); }
        part->orientations[type] = orientations;
      }
      if (inv_orient_map != NULL) {
        FmsIntMapCopy(orient_type, ent_orients, FMS_ORIENTATION_INT_TYPE,
                      orientations + num_entities, inv_orient_map, num_ents);
      } else {
        FmsIntConvertCopy(orient_type, ent_orients, FMS_ORIENTATION_INT_TYPE,
                          orientations + num_entities, num_ents);
      }
    } else {
      // part->orientations[type] is already NULL since num_entities is 0
    }
  }
  // Update the part
  part->num_entities[type] = num_entities + num_ents;
  if (num_entities == 0) {
    part->entities_ids_type[type] = id_store_type;
  }
  // Update the component
  if (comp->dim == FMS_INVALID_DIM) {
    comp->dim = FmsEntityDim[type];
  }
  comp->num_main_entities += num_ents;
  return 0;
}

int FmsComponentAddPartSubEntities(FmsComponent comp, FmsInt part_id,
                                   FmsEntityType type, FmsIntType id_store_type,
                                   FmsIntType ent_id_type, const void *ents,
                                   FmsInt num_ents) {
  if (!comp) { E_RETURN(1); }
  if (part_id >= comp->num_parts) { E_RETURN(2); }
  if (type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(3); }
  if (comp->dim != FMS_INVALID_DIM && comp->dim <= FmsEntityDim[type]) {
    E_RETURN(4);
  }
  if (id_store_type >= FMS_NUM_INT_TYPES) { E_RETURN(5); }
  if (ent_id_type >= FMS_NUM_INT_TYPES) { E_RETURN(6); }
  struct _FmsPart_private *part = comp->parts + part_id;
  const FmsInt num_entities = part->num_entities[type];
  void *entities_ids = part->entities_ids[type];
  if (num_entities != 0) {
    if (entities_ids == NULL) { E_RETURN(7); }
    if (ents == NULL) { E_RETURN(8); }
    if (part->entities_ids_type[type] != id_store_type) { E_RETURN(9); }
  } else if (ents == NULL) {
    if (num_ents != part->domain->num_entities[type]) { E_RETURN(10); }
  }
  const size_t sizeof_ent_ids = FmsIntTypeSize[id_store_type];
  if (num_ents != 0) {
    const FmsInt nc = NextPow2(num_entities + num_ents);
    if (ents != NULL) {
      if (num_entities <= nc/2) {
        entities_ids = realloc(entities_ids, nc*sizeof_ent_ids);
        if (entities_ids == NULL) { E_RETURN(11); }
        part->entities_ids[type] = entities_ids;
      }
      FmsIntConvertCopy(ent_id_type, ents, id_store_type,
                        (char*)entities_ids + num_entities*sizeof_ent_ids,
                        num_ents);
    } else {
      // part->entities_ids[type] is already NULL since num_entities is 0
    }
  }
  // Update the part
  part->num_entities[type] = num_entities + num_ents;
  if (num_entities == 0) {
    part->entities_ids_type[type] = id_store_type;
  }
  // Nothing to update in the component
  return 0;
}

int FmsComponentAddRelation(FmsComponent comp, FmsInt other_comp_id) {
  if (!comp) { E_RETURN(1); }
  FmsInt num_relations = comp->num_relations;
  FmsInt *relations = comp->relations;
  FmsInt nc = NextPow2(num_relations + 1);
  if (num_relations <= nc/2) {
    relations = realloc(relations, nc*sizeof(FmsInt));
    if (relations == NULL) { E_RETURN(2); }
    comp->relations = relations;
  }
  relations[num_relations] = other_comp_id;
  comp->num_relations = num_relations + 1;
  return 0;
}

int FmsComponentSetCoordinates(FmsComponent comp, FmsField coords) {
  if (!comp) { E_RETURN(1); }
  if (coords->fd->component != comp) { E_RETURN(2); }
  comp->coordinates = coords;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsTag functions: construction */
/* -------------------------------------------------------------------------- */

int FmsTagSetComponent(FmsTag tag, FmsComponent comp) {
  if (!tag) { E_RETURN(1); }
  if (!comp) { E_RETURN(2); }
  tag->comp = comp;
  return 0;
}

int FmsTagSet(FmsTag tag, FmsIntType stored_tag_type, FmsIntType input_tag_type,
              const void *ent_tags, FmsInt num_ents) {
  if (!tag) { E_RETURN(1); }
  if (stored_tag_type >= FMS_NUM_INT_TYPES) { E_RETURN(2); }
  if (tag->tag_type != FMS_NUM_INT_TYPES && tag->tag_type != stored_tag_type) {
    E_RETURN(3);
  }
  if (input_tag_type >= FMS_NUM_INT_TYPES) { E_RETURN(4); }
  FmsComponent comp = tag->comp;
  if (!comp) { E_RETURN(5); }
  if (comp->num_main_entities != num_ents) { E_RETURN(6); }
  if (num_ents != 0 && ent_tags == NULL) { E_RETURN(7); }
  free(tag->tags);
  tag->tags = NULL;
  if (num_ents != 0) {
    tag->tags = malloc(num_ents*FmsIntTypeSize[stored_tag_type]);
    if (tag->tags == NULL) { E_RETURN(8); }
    FmsIntConvertCopy(input_tag_type, ent_tags, stored_tag_type, tag->tags,
                      num_ents);
  }
  tag->tag_type = stored_tag_type;
  return 0;
}

int FmsTagAllocate(FmsTag tag, FmsIntType stored_tag_type, void **ent_tags,
                   FmsInt *num_ents) {
  if (!tag) { E_RETURN(1); }
  if (stored_tag_type >= FMS_NUM_INT_TYPES) { E_RETURN(2); }
  if (tag->tag_type != FMS_NUM_INT_TYPES && tag->tag_type != stored_tag_type) {
    E_RETURN(3);
  }
  if (ent_tags == NULL) { E_RETURN(4); }
  FmsComponent comp = tag->comp;
  if (!comp) { E_RETURN(5); }
  free(tag->tags);
  tag->tags = NULL;
  FmsInt num_entities = comp->num_main_entities;
  if (num_entities != 0) {
    tag->tags = malloc(num_entities*FmsIntTypeSize[stored_tag_type]);
    if (tag->tags == NULL) { E_RETURN(6); }
  }
  tag->tag_type = stored_tag_type;
  *ent_tags = tag->tags;
  if (num_ents != NULL) { *num_ents = num_entities; }
  return 0;
}

int FmsTagAddDescriptions(FmsTag tag, FmsIntType tag_type, const void *tags,
                          const char *const *tag_descr, FmsInt num_tags) {
  if (!tag) { E_RETURN(1); }
  if (tag_type >= FMS_NUM_INT_TYPES) { E_RETURN(2); }
  if (tag->tag_type == FMS_NUM_INT_TYPES) { E_RETURN(3); }
  if (num_tags == 0) { return 0; }
  if (tags == NULL) { E_RETURN(4); }
  if (tag_descr == NULL) { E_RETURN(5); }
  const FmsInt num_tag_descriptions = tag->num_tag_descriptions;
  char *described_tags = tag->described_tags;
  char **tag_descriptions = tag->tag_descriptions;
  const size_t sizeof_tag_type = FmsIntTypeSize[tag->tag_type];
  const FmsInt nc = NextPow2(num_tag_descriptions + num_tags);
  if (num_tag_descriptions <= nc/2) {
    described_tags = realloc(described_tags, nc*sizeof_tag_type);
    if (described_tags == NULL) { E_RETURN(6); }
    tag->described_tags = described_tags;
    tag_descriptions = realloc(tag_descriptions, nc*sizeof(char*));
    if (tag_descriptions == NULL) { E_RETURN(6); }
    tag->tag_descriptions = tag_descriptions;
  }
  // Copy the new described tags
  FmsIntConvertCopy(tag_type, tags, tag->tag_type,
                    described_tags + num_tag_descriptions*sizeof_tag_type,
                    num_tags);
  // Copy the new description strings
  for (FmsInt i = 0; i < num_tags; i++) {
    if (FmsCopyString(tag_descr[i],
                      &tag_descriptions[num_tag_descriptions + i])) {
      while (i != 0) {
        free(tag_descriptions[num_tag_descriptions + (--i)]);
      }
      E_RETURN(6);
    }
  }
  // Update tag
  tag->num_tag_descriptions = num_tag_descriptions + num_tags;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsDataCollection functions: construction */
/* -------------------------------------------------------------------------- */

int FmsDataCollectionCreate(FmsMesh mesh, const char *dc_name,
                            FmsDataCollection *dc) {
  if (!mesh) { E_RETURN(1); }
  if (!dc) { E_RETURN(2); }
  FmsDataCollection dcoll;
  dcoll = calloc(1, sizeof(*dcoll));
  if (dcoll == NULL) { E_RETURN(3); }
  // On error, call free(dcoll)
  // Init dcoll
  if (FmsCopyString(dc_name, &dcoll->name)) { free(dcoll); E_RETURN(3); }
  dcoll->mesh = mesh;
  *dc = dcoll;
  return 0;
}

int FmsDataCollectionDestroy(FmsDataCollection *dc) {
  int err = 0;
  if (!dc) { E_RETURN(1); }
  FmsDataCollection dcoll = *dc;
  if (dcoll != NULL) {
    // Delete fields
    FmsInt num_fields = dcoll->num_fields;
    FmsField *fields = dcoll->fields;
    for (FmsInt i = 0; i < num_fields; i++) {
      FmsField field = fields[i];
      free(field->data);
      err = UpdateError(err, FmsMetaDataDestroy(&field->mdata));
      free(field->name);
      free(field);
    }
    free(fields);

    // Delete field-descriptors
    FmsInt num_fds = dcoll->num_fds;
    FmsFieldDescriptor *fds = dcoll->fds;
    for (FmsInt i = 0; i < num_fds; i++) {
      FmsFieldDescriptor fd = fds[i];
      switch (fd->descr_type) {
      case FMS_FIXED_ORDER: free(fd->descriptor.fixed_order); break;
      default: err = UpdateError(err, 2);
      }
      free(fd->name);
      free(fd);
    }
    free(fds);

    // Delete mesh, meta-data, and name
    err = UpdateError(err, FmsMeshDestroy(&dcoll->mesh));
    err = UpdateError(err, FmsMetaDataDestroy(&dcoll->mdata));
    free(dcoll->name);

    // Delete dcoll
    free(dcoll);
    *dc = NULL;
  }
  E_RETURN(err);
}

int FmsDataCollectionAddFieldDescriptor(FmsDataCollection dc,
                                        const char *fd_name,
                                        FmsFieldDescriptor *fd) {
  if (!dc) { E_RETURN(1); }
  if (!fd) { E_RETURN(2); }
  const FmsInt num_fds = dc->num_fds;
  FmsFieldDescriptor *fds = dc->fds;
  // Reallocate dc->fds, if necessary
  const FmsInt nc = NextPow2(num_fds + 1);
  if (num_fds <= nc/2) {
    fds = realloc(fds, nc*sizeof(*fds));
    if (fds == NULL) { E_RETURN(3); }
    dc->fds = fds;
  }
  FmsFieldDescriptor new_fd;
  new_fd = calloc(1, sizeof(*new_fd));
  if (new_fd == NULL) { E_RETURN(3); }
  // On error, call free(new_fd)
  if (FmsCopyString(fd_name, &new_fd->name)) { free(new_fd); E_RETURN(3); }
  // Update the field descriptor
  dc->num_fds = num_fds + 1;
  dc->fds[num_fds] = new_fd;
  *fd = new_fd;
  return 0;
}

int FmsDataCollectionAddField(FmsDataCollection dc, const char *field_name,
                              FmsField *field) {
  if (!dc) { E_RETURN(1); }
  if (!field) { E_RETURN(2); }
  const FmsInt num_fields = dc->num_fields;
  FmsField *fields = dc->fields;
  // Reallocate dc->fields, if necessary
  const FmsInt nc = NextPow2(num_fields + 1);
  if (num_fields <= nc/2) {
    fields = realloc(fields, nc*sizeof(*fields));
    if (fields == NULL) { E_RETURN(3); }
    dc->fields = fields;
  }
  FmsField new_field;
  new_field = calloc(1, sizeof(*new_field));
  if (new_field == NULL) { E_RETURN(3); }
  // On error, call free(new_field)
  if (FmsCopyString(field_name, &new_field->name)) {
    free(new_field); E_RETURN(3);
  }
  // Update the field
  dc->num_fields = num_fields + 1;
  dc->fields[num_fields] = new_field;
  *field = new_field;
  return 0;
}

int FmsDataCollectionAttachMetaData(FmsDataCollection dc, FmsMetaData *mdata) {
  if (!dc) { E_RETURN(1); }
  if (!mdata) { E_RETURN(2); }
  FmsMetaData md = dc->mdata;
  if (md == NULL) {
    md = calloc(1, sizeof(*md));
    if (md == NULL) { E_RETURN(3); }
    dc->mdata = md;
  } else {
    int err = FmsMetaDataClear(md);
    if (err) { E_RETURN(err); }
  }
  *mdata = md;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsFieldDescriptor functions: construction */
/* -------------------------------------------------------------------------- */

int FmsFieldDescriptorSetComponent(FmsFieldDescriptor fd, FmsComponent comp) {
  if (!fd) { E_RETURN(1); }
  fd->component = comp;
  return 0;
}

int FmsFieldDescriptorSetFixedOrder(FmsFieldDescriptor fd,
                                    FmsFieldType field_type,
                                    FmsBasisType basis_type, FmsInt order) {
  if (!fd) { E_RETURN(1); }
  FmsComponent comp = fd->component;
  if (!comp) { E_RETURN(2); }
  if (fd->descriptor.any != NULL) { E_RETURN(3); }
  FmsInt num_dofs = 0;
  // Count the number of dofs - based on component, field_type and order
  FmsInt ent_dofs[FMS_NUM_ENTITY_TYPES];
  const FmsInt dim = comp->dim;
  const FmsInt p = order;
  switch (field_type) {
  case FMS_CONTINUOUS: {
    const FmsInt pm1 = p-1, pm2 = p-2, pm3 = p-3;
    if (p <= 0) { E_RETURN(4); }
    ent_dofs[FMS_VERTEX] = 1;
    ent_dofs[FMS_EDGE] = pm1;
    ent_dofs[FMS_TRIANGLE] = (pm1*pm2)/2;
    ent_dofs[FMS_QUADRILATERAL] = pm1*pm1;
    ent_dofs[FMS_TETRAHEDRON] = (ent_dofs[FMS_TRIANGLE]*pm3)/3;
    ent_dofs[FMS_HEXAHEDRON] = ent_dofs[FMS_QUADRILATERAL]*pm1;
    ent_dofs[FMS_WEDGE] = ent_dofs[FMS_TRIANGLE]*pm1;
    ent_dofs[FMS_PYRAMID] = (ent_dofs[FMS_TRIANGLE]*(p+pm3))/3;
    break;
  }
  case FMS_DISCONTINUOUS:
  case FMS_DISCONTINUOUS_WEIGHTED: {
    const FmsInt pp1 = p+1, pp2 = p+2, pp3 = p+3;
    // if (p < 0) { E_RETURN(4); } // If we make FmsInt a signed type
    ent_dofs[FMS_VERTEX] = 1;
    ent_dofs[FMS_EDGE] = pp1;
    ent_dofs[FMS_TRIANGLE] = (pp1*pp2)/2;
    ent_dofs[FMS_QUADRILATERAL] = pp1*pp1;
    ent_dofs[FMS_TETRAHEDRON] = (ent_dofs[FMS_TRIANGLE]*pp3)/3;
    ent_dofs[FMS_HEXAHEDRON] = ent_dofs[FMS_QUADRILATERAL]*pp1;
    ent_dofs[FMS_WEDGE] = ent_dofs[FMS_TRIANGLE]*pp1;
    ent_dofs[FMS_PYRAMID] = (ent_dofs[FMS_TRIANGLE]*(p+pp3))/3;
    for (FmsInt et = 0; et < FMS_NUM_ENTITY_TYPES; et++) {
      if (FmsEntityDim[et] < dim) {
        ent_dofs[et] = 0;
      }
    }
    break;
  }
  case FMS_HCURL: {
    FmsAbortNotImplemented();
    break;
  }
  case FMS_HDIV: {
    FmsAbortNotImplemented();
    break;
  }
  default: E_RETURN(5);
  }
  // Count the dofs part by part
  for (FmsInt i = 0; i < comp->num_parts; i++) {
    for (FmsInt et = 0; et < FMS_NUM_ENTITY_TYPES; et++) {
      num_dofs += comp->parts[i].num_entities[et] * ent_dofs[et];
    }
  }
  // Construct fixed-order descriptor
  struct _FmsFieldDescriptor_FixedOrder *fo;
  fo = malloc(sizeof(*fo));
  if (fo == NULL) { E_RETURN(6); }
  fo->field_type = field_type;
  fo->basis_type = basis_type;
  fo->order = order;
  // Update fd
  fd->descr_type = FMS_FIXED_ORDER;
  fd->descriptor.fixed_order = fo;
  fd->num_dofs = num_dofs;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsField functions: construction */
/* -------------------------------------------------------------------------- */

int FmsFieldSet(FmsField field, FmsFieldDescriptor fd, FmsInt num_vec_comp,
                FmsLayoutType layout_type, FmsScalarType data_type,
                const void *data) {
  if (!field) { E_RETURN(1); }
  if (field->fd != NULL) { E_RETURN(2); }
  if (!fd) { E_RETURN(3); }
  if (!fd->descriptor.any) { E_RETURN(4); }
  if (layout_type != FMS_BY_NODES &&
      layout_type != FMS_BY_VDIM) { E_RETURN(5); }
  if (data_type >= FMS_NUM_SCALAR_TYPES) { E_RETURN(6); }
  if (!data) { E_RETURN(7); }
  // Make a copy of the data array
  const FmsInt num_vdofs = fd->num_dofs * num_vec_comp;
  const size_t sizeof_data_type = FmsScalarTypeSize[data_type];
  void *data_copy;
  if (num_vdofs != 0) {
    data_copy = malloc(num_vdofs*sizeof_data_type);
    if (data_copy == NULL) { E_RETURN(8); }
    memcpy(data_copy, data, num_vdofs*sizeof_data_type);
  } else {
    data_copy = NULL;
  }
  // Update field
  field->fd = fd;
  field->num_vec_comp = num_vec_comp;
  field->layout = layout_type;
  field->scalar_type = data_type;
  field->data = data_copy;
  return 0;
}

int FmsFieldAttachMetaData(FmsField field, FmsMetaData *mdata) {
  if (!field) { E_RETURN(1); }
  if (!mdata) { E_RETURN(2); }
  FmsMetaData md = field->mdata;
  if (md == NULL) {
    md = calloc(1, sizeof(*md));
    if (md == NULL) { E_RETURN(3); }
    field->mdata = md;
  } else {
    int err = FmsMetaDataClear(md);
    if (err) { E_RETURN(err); }
  }
  *mdata = md;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsMetaData functions: construction */
/* -------------------------------------------------------------------------- */

int FmsMetaDataSetIntegers(FmsMetaData mdata, const char *mdata_name,
                           FmsIntType int_type, FmsInt size, void **data) {
  if (!mdata) { E_RETURN(1); }
  if (int_type >= FMS_NUM_INT_TYPES) { E_RETURN(2); }
  if (size && !data) { E_RETURN(3); }
  int err = FmsMetaDataClear(mdata);
  if (err) { E_RETURN(err); }
  const size_t sizeof_int_type = FmsIntTypeSize[int_type];
  void *md_data = size ? malloc(size*sizeof_int_type) : NULL;
  if (size && md_data == NULL) { E_RETURN(4); }
  // On error, call free(md_data)
  // Update mdata
  if (FmsCopyString(mdata_name, &mdata->name)) { free(md_data); E_RETURN(4); }
  mdata->md_type = FMS_INTEGER;
  mdata->sub_type.int_type = int_type;
  mdata->num_entries = size;
  mdata->data = md_data;
  *data = md_data;
  return 0;
}

int FmsMetaDataSetScalars(FmsMetaData mdata, const char *mdata_name,
                          FmsScalarType scal_type, FmsInt size, void **data) {
  if (!mdata) { E_RETURN(1); }
  if (scal_type >= FMS_NUM_SCALAR_TYPES) { E_RETURN(2); }
  if (size && !data) { E_RETURN(3); }
  int err = FmsMetaDataClear(mdata);
  if (err) { E_RETURN(err); }
  const size_t sizeof_scal_type = FmsScalarTypeSize[scal_type];
  void *md_data = size ? malloc(size*sizeof_scal_type) : NULL;
  if (size && md_data == NULL) { E_RETURN(4); }
  // On error, call free(md_data)
  // Update mdata
  if (FmsCopyString(mdata_name, &mdata->name)) { free(md_data); E_RETURN(4); }
  mdata->md_type = FMS_SCALAR;
  mdata->sub_type.scalar_type = scal_type;
  mdata->num_entries = size;
  mdata->data = md_data;
  *data = md_data;
  return 0;
}

int FmsMetaDataSetString(FmsMetaData mdata, const char *mdata_name,
                         const char *c_string) {
  if (!mdata) { E_RETURN(1); }
  int err = FmsMetaDataClear(mdata);
  if (err) { E_RETURN(err); }
  char *md_data;
  if (FmsCopyString(c_string, &md_data)) { E_RETURN(2); }
  // On error, call free(md_data)
  // Update mdata
  if (FmsCopyString(mdata_name, &mdata->name)) { free(md_data); E_RETURN(2); }
  mdata->md_type = FMS_STRING;
  mdata->num_entries = md_data ? strlen(md_data) : 0;
  mdata->data = md_data;
  return 0;
}

int FmsMetaDataSetMetaData(FmsMetaData mdata, const char *mdata_name,
                           FmsInt size, FmsMetaData **data) {
  if (!mdata) { E_RETURN(1); }
  int err = FmsMetaDataClear(mdata);
  if (err) { E_RETURN(err); }
  char *name_copy;
  if (FmsCopyString(mdata_name, &name_copy)) { E_RETURN(2); }
  // On error, call free(name_copy)
  FmsMetaData *md_data = size ? malloc(size*sizeof(*md_data)) : NULL;
  if (size && md_data == NULL) { free(name_copy); E_RETURN(2); }
  // On error, call free(md_data)
  for (FmsInt i = 0; i < size; i++) {
    FmsMetaData md;
    md = calloc(1, sizeof(*md));
    if (md == NULL) {
      while (i != 0) { free(md_data[--i]); }
      free(md_data);
      free(name_copy);
      E_RETURN(2);
    }
    md_data[i] = md;
  }
  // Update mdata
  mdata->name = name_copy;
  mdata->md_type = FMS_META_DATA;
  mdata->num_entries = size;
  mdata->data = md_data;
  *data = md_data;
  return 0;
}

int FmsMetaDataClear(FmsMetaData mdata) {
  int err = 0;
  if (!mdata) { E_RETURN(1); }
  if (mdata->data) {
    switch (mdata->md_type) {
    case FMS_INTEGER:
    case FMS_SCALAR:
    case FMS_STRING: break;
    case FMS_META_DATA: {
      FmsMetaData *mds = (FmsMetaData*)(mdata->data);
      for (FmsInt i = mdata->num_entries; i != 0; ) {
        err = UpdateError(err, FmsMetaDataDestroy(&mds[--i]));
      }
      break;
    }
    default: err = UpdateError(err, 2);
    }
    free(mdata->data);
    mdata->data = NULL;
  }
  free(mdata->name);
  mdata->name = NULL;
  mdata->num_entries = 0;
  E_RETURN(err);
}

int FmsMetaDataDestroy(FmsMetaData *mdata) {
  int err = 0;
  if (!mdata) { E_RETURN(1); }
  FmsMetaData md = *mdata;
  if (md) {
    err = UpdateError(err, FmsMetaDataClear(md));
    free(md);
    *mdata = NULL;
  }
  E_RETURN(err);
}


/* -------------------------------------------------------------------------- */
/* FmsMesh functions: query interface */
/* -------------------------------------------------------------------------- */

int FmsMeshGetPartitionId(FmsMesh mesh, FmsInt *partition_id,
                          FmsInt *num_partitions) {
  if (!mesh) { E_RETURN(1); }
  if (partition_id) { *partition_id = mesh->partition_id; }
  if (num_partitions) { *num_partitions = mesh->num_partitions; }
  return 0;
}

int FmsMeshGetNumDomainNames(FmsMesh mesh, FmsInt *num_domain_names) {
  if (!mesh) { E_RETURN(1); }
  if (!num_domain_names) { E_RETURN(2); }
  *num_domain_names = mesh->num_domain_names;
  return 0;
}

int FmsMeshGetDomains(FmsMesh mesh, FmsInt domain_name_id,
                      const char **domain_name, FmsInt *num_domains,
                      FmsDomain **domains) {
  if (!mesh) { E_RETURN(1); }
  if (domain_name_id >= mesh->num_domain_names) { E_RETURN(2); }
  if (domain_name) { *domain_name = mesh->domain_names[domain_name_id]; }
  if (num_domains) {
    *num_domains = mesh->domain_offsets[domain_name_id+1] -
                   mesh->domain_offsets[domain_name_id];
  }
  if (domains) {
    *domains = mesh->domains + mesh->domain_offsets[domain_name_id];
  }
  return 0;
}

int FmsMeshGetDomainsByName(FmsMesh mesh, const char *domain_name,
                            FmsInt *num_domains, FmsDomain **domains) {
  if (!mesh) { E_RETURN(1); }
  if (!domain_name) { E_RETURN(2); }
  const FmsInt num_domain_names = mesh->num_domain_names;
  char **domain_names = mesh->domain_names;
  FmsInt domain_name_id = 0;
  for ( ; 1; domain_name_id++) {
    if (domain_name_id == num_domain_names) { E_RETURN(4); }
    if (strcmp(domain_name, domain_names[domain_name_id]) == 0) { break; }
  }
  if (num_domains) {
    *num_domains = mesh->domain_offsets[domain_name_id+1] -
                   mesh->domain_offsets[domain_name_id];
  }
  if (domains) {
    *domains = mesh->domains + mesh->domain_offsets[domain_name_id];
  }
  return 0;
}

int FmsMeshGetNumComponents(FmsMesh mesh, FmsInt *num_comp) {
  if (!mesh) { E_RETURN(1); }
  if (!num_comp) { E_RETURN(2); }
  *num_comp = mesh->num_comps;
  return 0;
}

int FmsMeshGetComponent(FmsMesh mesh, FmsInt comp_id, FmsComponent *comp) {
  if (!mesh) { E_RETURN(1); }
  if (comp_id >= mesh->num_comps) { E_RETURN(2); }
  if (!comp) { E_RETURN(3); }
  *comp = mesh->comps[comp_id];
  return 0;
}

int FmsMeshGetNumTags(FmsMesh mesh, FmsInt *num_tags) {
  if (!mesh) { E_RETURN(1); }
  if (!num_tags) { E_RETURN(2); }
  *num_tags = mesh->num_tags;
  return 0;
}

int FmsMeshGetTag(FmsMesh mesh, FmsInt tag_id, FmsTag *tag) {
  if (!mesh) { E_RETURN(1); }
  if (tag_id >= mesh->num_tags) { E_RETURN(2); }
  if (!tag) { E_RETURN(3); }
  *tag = mesh->tags[tag_id];
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsDomain functions: query interface */
/* -------------------------------------------------------------------------- */

int FmsDomainGetName(FmsDomain domain, const char **domain_name,
                     FmsInt *domain_id) {
  if (!domain) { E_RETURN(1); }
  if (domain_name) { *domain_name = domain->name; }
  if (domain_id) { *domain_id = domain->id; }
  return 0;
}

int FmsDomainGetDimension(FmsDomain domain, FmsInt *dim) {
  if (!domain) { E_RETURN(1); }
  if (!dim) { E_RETURN(2); }
  *dim = domain->dim;
  return 0;
}

int FmsDomainGetNumVertices(FmsDomain domain, FmsInt *num_verts) {
  if (!domain) { E_RETURN(1); }
  if (!num_verts) { E_RETURN(2); }
  *num_verts = domain->num_entities[FMS_VERTEX];
  return 0;
}

int FmsDomainGetNumEntities(FmsDomain domain, FmsEntityType type,
                            FmsInt *num_ents) {
  if (!domain) { E_RETURN(1); }
  if (type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  if (!num_ents) { E_RETURN(3); }
  *num_ents = domain->num_entities[type];
  return 0;
}

int FmsDomainGetAllEntities(FmsDomain domain, FmsEntityType type,
                            FmsIntType *ent_id_type, const void **ents,
                            FmsInt *num_ents) {
  if (!domain) { E_RETURN(1); }
  if (type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  if (ent_id_type) { *ent_id_type = domain->side_ids_type[type]; }
  if (ents) { *ents = domain->entities[type]; }
  if (num_ents) { *num_ents = domain->num_entities[type]; }
  return 0;
}

int FmsDomainGetEntities(FmsDomain domain, FmsEntityType type,
                         FmsEntityReordering reordering, FmsIntType ent_id_type,
                         FmsInt first_ent, void *ents, FmsInt num_ents) {
  if (!domain) { E_RETURN(1); }
  if (type == FMS_VERTEX || type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  if (ent_id_type >= FMS_NUM_INT_TYPES) { E_RETURN(3); }
  if (num_ents == 0) { return 0; }
  if (first_ent >= domain->num_entities[type]) { E_RETURN(4); }
  if (!ents) { E_RETURN(5); }
  if (first_ent + num_ents > domain->num_entities[type]) { E_RETURN(6); }
  void *src_ents = domain->entities[type];
  if (!src_ents) { E_RETURN(7); }
  const FmsIntType side_ids_type = domain->side_ids_type[type];
  const size_t sizeof_side_ids = FmsIntTypeSize[side_ids_type];
  const FmsInt num_sides = FmsEntityNumSides[type];
  src_ents = (char*)src_ents + first_ent*num_sides*sizeof_side_ids;
  if (reordering == NULL || reordering[type] == NULL) {
    FmsIntConvertCopy(side_ids_type, src_ents, ent_id_type, ents,
                      num_ents*num_sides);
  } else {
    const int *perm = reordering[type];
    int inv_perm[num_sides]; // variable length array
    for (FmsInt i = 0; i < num_sides; i++) {
      inv_perm[perm[i]] = i;
    }
    FmsIntPermuteConvertCopy(side_ids_type, src_ents, ent_id_type, ents,
                             inv_perm, num_sides, num_ents*num_sides);
  }
  return 0;
}

int FmsDomainGetEntitiesVerts(FmsDomain domain, FmsEntityType type,
                              FmsEntityReordering vert_reordering,
                              FmsIntType ent_id_type, FmsInt first_ent,
                              void *ents_verts, FmsInt num_ents) {
  if (!domain) { E_RETURN(1); }
  if (type == FMS_VERTEX || type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  if (ent_id_type >= FMS_NUM_INT_TYPES) { E_RETURN(3); }
  if (num_ents == 0) { return 0; }
  if (first_ent >= domain->num_entities[type]) { E_RETURN(4); }
  if (!ents_verts) { E_RETURN(5); }
  if (first_ent + num_ents > domain->num_entities[type]) { E_RETURN(6); }
  const char *src_ents = domain->entities[type];
  if (!src_ents) { E_RETURN(7); }
  const FmsIntType idx_type = domain->side_ids_type[FMS_EDGE];
  if (idx_type != domain->side_ids_type[type]) { E_RETURN(8); }
  const FmsOrientation *side_ori = domain->orientations[type];
  if (type != FMS_EDGE && !side_ori) { E_RETURN(9); }
  const int *perm = vert_reordering ? vert_reordering[type] : NULL;
  const size_t sizeof_idx_type = FmsIntTypeSize[idx_type];
  const FmsInt num_verts = FmsEntityNumVerts[type];
  const FmsInt num_sides = FmsEntityNumSides[type];
  src_ents += first_ent*num_sides*sizeof_idx_type;
  if (type != FMS_EDGE) { side_ori += first_ent*num_sides; }
  switch (type) {
  case FMS_EDGE:
    if (!perm) {
      FmsIntConvertCopy(idx_type, src_ents, ent_id_type, ents_verts,
                        num_ents*num_verts);
    } else {
      int inv_perm[2];
      for (FmsInt i = 0; i < 2; i++) { inv_perm[perm[i]] = i; }
      FmsIntPermuteConvertCopy(idx_type, src_ents, ent_id_type, ents_verts,
                               inv_perm, num_verts, num_ents*num_verts);
    }
    break;
  case FMS_TRIANGLE: {
    const void *all_edges = domain->entities[FMS_EDGE];
    int inv_perm[3];
    if (perm) {
      for (FmsInt i = 0; i < 3; i++) { inv_perm[perm[i]] = i; }
    }
    const int bsize = 1024;
    FmsInt tri_vert[4*bsize]; // the vertices of the first two edges
    FmsInt tri_vert_x[3*bsize]; // vertices before applying the permutation
    for (FmsInt off = 0; off < num_ents; off += bsize) {
      const FmsInt sz = FmsIntMin(bsize, num_ents-off);
      // src_ents -> tri_vert
      FmsEntitiesToSides(idx_type, num_sides, 0, 2, src_ents, 2,
                         all_edges, 4, 0, tri_vert, sz);
      // tri_vert + side_ori -> tri_vert_x
      for (FmsInt i = 0; i < sz; i++) {
        const FmsInt *tv = tri_vert+4*i;
        const FmsOrientation *eo = side_ori+3*i;
        FmsInt *verts = tri_vert_x+3*i;
        if (eo[0] == 0) {
          verts[0] = tv[0];
          verts[1] = tv[1];
        } else {
          verts[0] = tv[1];
          verts[1] = tv[0];
        }
        verts[2] = (eo[1] == 0) ? tv[3] : tv[2];
      }
      // tri_vert_x -> ents_verts, i.e. convert + vertex permutation
      if (!perm) {
        FmsIntConvertCopy(FMS_INT_TYPE, tri_vert_x, ent_id_type, ents_verts,
                          sz*num_verts);
      } else {
        FmsIntPermuteConvertCopy(FMS_INT_TYPE, tri_vert_x, ent_id_type,
                                 ents_verts, inv_perm, num_verts, sz*num_sides);
      }
      src_ents += sz*num_sides*sizeof_idx_type;
      side_ori += sz*num_sides;
      ents_verts = (char*)ents_verts + sz*num_verts*FmsIntTypeSize[ent_id_type];
    }
    break;
  }
  case FMS_QUADRILATERAL: {
    const void *all_edges = domain->entities[FMS_EDGE];
    int inv_perm[4];
    if (perm) {
      for (FmsInt i = 0; i < 4; i++) { inv_perm[perm[i]] = i; }
    }
    const int bsize = 1024;
    FmsInt quad_vert[4*bsize]; // the vertices of edges 0 and 2
    FmsInt quad_vert_x[4*bsize]; // vertices before applying the permutation
    for (FmsInt off = 0; off < num_ents; off += bsize) {
      const FmsInt sz = FmsIntMin(bsize, num_ents-off);
      // src_ents -> quad_vert
      FmsEntitiesToSides(idx_type, num_sides, 0, 1, src_ents, 2,
                         all_edges, 4, 0, quad_vert, sz);
      FmsEntitiesToSides(idx_type, num_sides, 2, 1, src_ents, 2,
                         all_edges, 4, 2, quad_vert, sz);
      // quad_vert + side_ori -> quad_vert_x
      for (FmsInt i = 0; i < sz; i++) {
        const FmsInt *qv = quad_vert+4*i;
        const FmsOrientation *eo = side_ori+4*i;
        FmsInt *verts = quad_vert_x+4*i;
        if (eo[0] == 0) {
          verts[0] = qv[0];
          verts[1] = qv[1];
        } else {
          verts[0] = qv[1];
          verts[1] = qv[0];
        }
        if (eo[2] == 0) {
          verts[2] = qv[2];
          verts[3] = qv[3];
        } else {
          verts[2] = qv[3];
          verts[3] = qv[2];
        }
      }
      // quad_vert_x -> ents_verts, i.e. convert + vertex permutation
      if (!perm) {
        FmsIntConvertCopy(FMS_INT_TYPE, quad_vert_x, ent_id_type, ents_verts,
                          sz*num_verts);
      } else {
        FmsIntPermuteConvertCopy(FMS_INT_TYPE, quad_vert_x, ent_id_type,
                                 ents_verts, inv_perm, num_verts, sz*num_sides);
      }
      src_ents += sz*num_sides*sizeof_idx_type;
      side_ori += sz*num_sides;
      ents_verts = (char*)ents_verts + sz*num_verts*FmsIntTypeSize[ent_id_type];
    }
    break;
  }
  case FMS_TETRAHEDRON: {
    if (idx_type != domain->side_ids_type[FMS_TRIANGLE]) { E_RETURN(10); }
    const void *all_edges = domain->entities[FMS_EDGE];
    const void *all_tris = domain->entities[FMS_TRIANGLE];
    const FmsOrientation *tri_side_ori = domain->orientations[FMS_TRIANGLE];
    if (!tri_side_ori) { E_RETURN(11); }
    int inv_perm[4];
    if (perm) {
      for (FmsInt i = 0; i < 4; i++) { inv_perm[perm[i]] = i; }
    }
    const int bsize = 1024;
    FmsInt tet_face[4*bsize]; // the 4 faces of the tets
    FmsInt tet_edge[6*bsize]; // the edges of faces B and C
    FmsInt tet_vert_x[4*bsize]; // vertices before applying the permutation
    for (FmsInt off = 0; off < num_ents; off += bsize) {
      const FmsInt sz = FmsIntMin(bsize, num_ents-off);
      // src_ents -> tet_face
      FmsIntConvertCopy(idx_type, src_ents, FMS_INT_TYPE, tet_face, sz*4);
      // src_ents -> tet_edge
      FmsEntitiesToSides(idx_type, num_sides, 1, 2, src_ents, 3,
                         all_tris, 6, 0, tet_edge, sz);
      // tet_edge -> tet_edge[mod]
      for (FmsInt i = 0; i < sz; i++) {
        const FmsInt *tf = tet_face+4*i;
        const FmsInt *te = tet_edge+6*i;
        const FmsOrientation *so = side_ori+4*i;
        FmsInt eo[6]; // side (edge) orientations for faces B and C
        FmsInt te_x[6], eo_x[6];
        for (int j = 0; j < 2; j++) {
          for (int e = 0; e < 3; e++) {
            eo[3*j+e] = tri_side_ori[3*tf[j+1]+e];
          }
        }
        if (so[1] == FMS_ORIENTATION_UNKNOWN ||
            so[2] == FMS_ORIENTATION_UNKNOWN) {
          E_RETURN(12);
        }
        FmsPermuteTriBdr(so[1], te+0, te_x+0);
        FmsPermuteTriBdr(so[1], eo+0, eo_x+0);
        if (so[1]%2) { for (int j = 0; j < 3; j++) { eo_x[j] = 1-eo_x[j]; } }
        FmsPermuteTriBdr(so[2], te+3, te_x+3);
        FmsPermuteTriBdr(so[2], eo+3, eo_x+3);
        if (so[2]%2) { for (int j = 3; j < 6; j++) { eo_x[j] = 1-eo_x[j]; } }
#ifndef NDEBUG
        if (te_x[2] != te_x[4]) { FmsInternalError(); }
        if (eo_x[2] != 1-eo_x[4]) { FmsInternalError(); }
#endif
        tet_edge[2*i+0] = te_x[0];
        tet_edge[2*i+1] = te_x[5];
        tet_face[2*i+0] = eo_x[0];
        tet_face[2*i+1] = eo_x[5];
      }
      // tet_edge[mod] -> tet_vert_x[mod]
      FmsIntConvertCopy(FMS_INT_TYPE, tet_edge,
                        idx_type, &tet_edge[2*bsize], 2*sz);
      FmsEntitiesToSides(idx_type, 2*sz, 0, 2*sz, &tet_edge[2*bsize],
                         2, all_edges, 4*sz, 0, tet_vert_x, 1);
      // tet_vert_x[mod] -> tet_vert_x
      for (FmsInt i = 0; i < sz; i++) {
        FmsInt *tv = tet_vert_x+4*i, z;
        if (tet_face[2*i+0] != 0) {
          z = tv[0];
          tv[0] = tv[1];
          tv[1] = z;
        }
        if (tet_face[2*i+1] == 0) {
          z = tv[2];
          tv[2] = tv[3];
          tv[3] = z;
        }
      }
      // tet_vert_x -> ents_verts, i.e. convert + vertex permutation
      if (!perm) {
        FmsIntConvertCopy(FMS_INT_TYPE, tet_vert_x, ent_id_type, ents_verts,
                          sz*num_verts);
      } else {
        FmsIntPermuteConvertCopy(FMS_INT_TYPE, tet_vert_x, ent_id_type,
                                 ents_verts, inv_perm, num_verts, sz*num_sides);
      }
      src_ents += sz*num_sides*sizeof_idx_type;
      side_ori += sz*num_sides;
      ents_verts = (char*)ents_verts + sz*num_verts*FmsIntTypeSize[ent_id_type];
    }
    break;
  }
  case FMS_HEXAHEDRON: {
    if (idx_type != domain->side_ids_type[FMS_QUADRILATERAL]) { E_RETURN(10); }
    const void *all_edges = domain->entities[FMS_EDGE];
    const void *all_quads = domain->entities[FMS_QUADRILATERAL];
    const FmsOrientation *quad_side_ori =
      domain->orientations[FMS_QUADRILATERAL];
    if (!quad_side_ori) { E_RETURN(11); }
    int inv_perm[8];
    if (perm) {
      for (FmsInt i = 0; i < 8; i++) { inv_perm[perm[i]] = i; }
    }
    const int bsize = 1024;
    FmsInt hex_face[6*bsize]; // the 4 faces of the tets
    FmsInt hex_edge[8*bsize]; // the edges of faces A and B
    FmsInt hex_vert_x[8*bsize]; // vertices before applying the permutation
    for (FmsInt off = 0; off < num_ents; off += bsize) {
      const FmsInt sz = FmsIntMin(bsize, num_ents-off);
      // src_ents -> hex_face
      FmsIntConvertCopy(idx_type, src_ents, FMS_INT_TYPE, hex_face, sz*6);
      // src_ents -> hex_edge
      FmsEntitiesToSides(idx_type, num_sides, 0, 2, src_ents, 4,
                         all_quads, 8, 0, hex_edge, sz);
      // hex_edge -> hex_edge[mod]
      for (FmsInt i = 0; i < sz; i++) {
        const FmsInt *hf = hex_face+6*i;
        const FmsInt *he = hex_edge+8*i;
        const FmsOrientation *so = side_ori+6*i;
        FmsInt eo[8]; // side (edge) orientations for faces A and B
        FmsInt he_x[8], eo_x[8];
        for (int j = 0; j < 2; j++) {
          for (int e = 0; e < 4; e++) {
            eo[4*j+e] = quad_side_ori[4*hf[j+1]+e];
          }
        }
        if (so[0] == FMS_ORIENTATION_UNKNOWN ||
            so[1] == FMS_ORIENTATION_UNKNOWN) {
          E_RETURN(12);
        }
        FmsPermuteQuadBdr(so[0], he+0, he_x+0);
        FmsPermuteQuadBdr(so[0], eo+0, eo_x+0);
        if (so[0]%2) { for (int j = 0; j < 4; j++) { eo_x[j] = 1-eo_x[j]; } }
        FmsPermuteQuadBdr(so[1], he+4, he_x+4);
        FmsPermuteQuadBdr(so[1], eo+4, eo_x+4);
        if (so[1]%2) { for (int j = 4; j < 8; j++) { eo_x[j] = 1-eo_x[j]; } }
        hex_edge[4*i+0] = he_x[0];
        hex_edge[4*i+1] = he_x[2];
        hex_edge[4*i+2] = he_x[4];
        hex_edge[4*i+3] = he_x[6];
        hex_face[4*i+0] = eo_x[0];
        hex_face[4*i+1] = eo_x[2];
        hex_face[4*i+2] = eo_x[4];
        hex_face[4*i+3] = eo_x[6];
      }
      // hex_edge[mod] -> hex_vert_x[mod]
      FmsIntConvertCopy(FMS_INT_TYPE, hex_edge,
                        idx_type, &hex_edge[4*bsize], 4*sz);
      FmsEntitiesToSides(idx_type, 4*sz, 0, 4*sz, &hex_edge[4*bsize],
                         2, all_edges, 8*sz, 0, hex_vert_x, 1);
      // hex_vert_x[mod] -> tet_vert_x
      for (FmsInt i = 0; i < sz; i++) {
        FmsInt *hv = hex_vert_x+8*i, z;
        if (hex_face[4*i+0] == 0) {
          z = hv[0];
          hv[0] = hv[1];
          hv[1] = z;
        }
        if (hex_face[4*i+1] == 0) {
          z = hv[2];
          hv[2] = hv[3];
          hv[3] = z;
        }
        if (hex_face[4*i+2] != 0) {
          z = hv[4];
          hv[4] = hv[5];
          hv[5] = z;
        }
        if (hex_face[4*i+3] != 0) {
          z = hv[6];
          hv[6] = hv[7];
          hv[7] = z;
        }
      }
      // hex_vert_x -> ents_verts, i.e. convert + vertex permutation
      if (!perm) {
        FmsIntConvertCopy(FMS_INT_TYPE, hex_vert_x, ent_id_type, ents_verts,
                          sz*num_verts);
      } else {
        FmsIntPermuteConvertCopy(FMS_INT_TYPE, hex_vert_x, ent_id_type,
                                 ents_verts, inv_perm, num_verts, sz*num_sides);
      }
      src_ents += sz*num_sides*sizeof_idx_type;
      side_ori += sz*num_sides;
      ents_verts = (char*)ents_verts + sz*num_verts*FmsIntTypeSize[ent_id_type];
    }
    break;
  }
  case FMS_WEDGE:
  case FMS_PYRAMID:
  default: FmsAbortNotImplemented();
  }
  return 0;
}

int FmsDomainGetAllOrientations(FmsDomain domain, FmsEntityType type,
                                const FmsOrientation **side_orients,
                                FmsInt *num_ents) {
  if (!domain) { E_RETURN(1); }
  if (type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(2); }
  if (side_orients) { *side_orients = domain->orientations[type]; }
  if (num_ents) { *num_ents = domain->num_entities[type]; }
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsComponent functions: query interface */
/* -------------------------------------------------------------------------- */

int FmsComponentGetName(FmsComponent comp, const char **comp_name) {
  if (!comp) { E_RETURN(1); }
  if (comp_name) { *comp_name = comp->name; }
  return 0;
}

int FmsComponentGetDimension(FmsComponent comp, FmsInt *dim) {
  if (!comp) { E_RETURN(1); }
  if (dim) { *dim = comp->dim; }
  return 0;
}

int FmsComponentGetNumParts(FmsComponent comp, FmsInt *num_parts) {
  if (!comp) { E_RETURN(1); }
  if (num_parts) { *num_parts = comp->num_parts; }
  return 0;
}

int FmsComponentGetPart(FmsComponent comp, FmsInt part_id, FmsEntityType type,
                        FmsDomain *domain, FmsIntType *ent_id_type,
                        const void **ents, const FmsOrientation **ent_orients,
                        FmsInt *num_ents) {
  if (!comp) { E_RETURN(1); }
  if (part_id >= comp->num_parts) { E_RETURN(2); }
  if (type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(3); }
  // if (FmsEntityDim[type] != comp->dim) { E_RETURN(4); }
  struct _FmsPart_private *part = &comp->parts[part_id];
  if (domain) { *domain = part->domain; }
  if (ent_id_type) { *ent_id_type = part->entities_ids_type[type]; }
  if (ents) { *ents = part->entities_ids[type]; }
  if (ent_orients) { *ent_orients = part->orientations[type]; }
  if (num_ents) { *num_ents = part->num_entities[type]; }
  return 0;
}

int FmsComponentGetPartSubEntities(FmsComponent comp, FmsInt part_id,
                                   FmsEntityType type, FmsIntType *ent_id_type,
                                   const void **ents, FmsInt *num_ents) {
  if (!comp) { E_RETURN(1); }
  if (part_id >= comp->num_parts) { E_RETURN(2); }
  if (type >= FMS_NUM_ENTITY_TYPES) { E_RETURN(3); }
  // if (FmsEntityDim[type] >= comp->dim) { E_RETURN(4); }
  struct _FmsPart_private *part = &comp->parts[part_id];
  if (ent_id_type) { *ent_id_type = part->entities_ids_type[type]; }
  if (ents) { *ents = part->entities_ids[type]; }
  if (num_ents) { *num_ents = part->num_entities[type]; }
  return 0;
}

int FmsComponentGetNumEntities(FmsComponent comp, FmsInt *num_ents) {
  if (!comp) { E_RETURN(1); }
  if (num_ents) { *num_ents = comp->num_main_entities; }
  return 0;
}

int FmsComponentGetRelations(FmsComponent comp, const FmsInt **rel_comps,
                             FmsInt *num_rel_comps) {
  if (!comp) { E_RETURN(1); }
  if (rel_comps) { *rel_comps = comp->relations; }
  if (num_rel_comps) { *num_rel_comps = comp->num_relations; }
  return 0;
}

int FmsComponentGetCoordinates(FmsComponent comp, FmsField *coords) {
  if (!comp) { E_RETURN(1); }
  if (coords) { *coords = comp->coordinates; }
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsTag functions: query interface */
/* -------------------------------------------------------------------------- */

int FmsTagGetName(FmsTag tag, const char **tag_name) {
  if (!tag) { E_RETURN(1); }
  if (tag_name) { *tag_name = tag->name; }
  return 0;
}

int FmsTagGetComponent(FmsTag tag, FmsComponent *comp) {
  if (!tag) { E_RETURN(1); }
  if (comp) { *comp = tag->comp; }
  return 0;
}

int FmsTagGet(FmsTag tag, FmsIntType *tag_type, const void **ent_tags,
              FmsInt *num_ents) {
  if (!tag) { E_RETURN(1); }
  if (tag_type) { *tag_type = tag->tag_type; }
  if (ent_tags) { *ent_tags = tag->tags; }
  if (num_ents) { *num_ents = tag->comp->num_main_entities; }
  return 0;
}

int FmsTagGetDescriptions(FmsTag tag, FmsIntType *tag_type, const void **tags,
                          const char *const **tag_descr, FmsInt *num_tags) {
  if (!tag) { E_RETURN(1); }
  if (tag_type) { *tag_type = tag->tag_type; }
  if (tags) { *tags = tag->described_tags; }
  if (tag_descr) {
    // without the explicit typecast, there is a warning; why?
    *tag_descr = (const char *const *)(tag->tag_descriptions);
  }
  if (num_tags) { *num_tags = tag->num_tag_descriptions; }
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsDataCollection functions: query interface */
/* -------------------------------------------------------------------------- */

int FmsDataCollectionGetName(FmsDataCollection dc, char **name) {
  if (!dc) { E_RETURN(1); }
  if (!name) { E_RETURN(2); }
  *name = dc->name;
  return 0;
}

int FmsDataCollectionGetMesh(FmsDataCollection dc, FmsMesh *mesh) {
  if (!dc) { E_RETURN(1); }
  if (!mesh) { E_RETURN(2); }
  *mesh = dc->mesh;
  return 0;
}

int FmsDataCollectionGetFieldDescriptors(FmsDataCollection dc,
    FmsFieldDescriptor **fds, FmsInt *num_fds) {
  if (!dc) { E_RETURN(1); }
  if (fds) { *fds = dc->fds; }
  if (num_fds) { *num_fds = dc->num_fds; }
  return 0;
}

int FmsDataCollectionGetFields(FmsDataCollection dc, FmsField **fields,
                               FmsInt *num_fields) {
  if (!dc) { E_RETURN(1); }
  if (fields) { *fields = dc->fields; }
  if (num_fields) { *num_fields = dc->num_fields; }
  return 0;
}

int FmsDataCollectionGetMetaData(FmsDataCollection dc, FmsMetaData *mdata) {
  if (!dc) { E_RETURN(1); }
  if (!mdata) { E_RETURN(2); }
  *mdata = dc->mdata;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsFieldDescriptor functions: query interface */
/* -------------------------------------------------------------------------- */

int FmsFieldDescriptorGetName(FmsFieldDescriptor fd, const char **fd_name) {
  if (!fd) { E_RETURN(1); }
  if (!fd_name) { E_RETURN(2); }
  *fd_name = fd->name;
  return 0;
}

int FmsFieldDescriptorGetComponent(FmsFieldDescriptor fd, FmsComponent *comp) {
  if (!fd) { E_RETURN(1); }
  if (!comp) { E_RETURN(2); }
  *comp = fd->component;
  return 0;
}

int FmsFieldDescriptorGetType(FmsFieldDescriptor fd,
                              FmsFieldDescriptorType *fd_type) {
  if (!fd) { E_RETURN(1); }
  if (!fd_type) { E_RETURN(2); }
  *fd_type = fd->descr_type;
  return 0;
}

int FmsFieldDescriptorGetFixedOrder(FmsFieldDescriptor fd,
                                    FmsFieldType *field_type,
                                    FmsBasisType *basis_type, FmsInt *order) {
  if (!fd) { E_RETURN(1); }
  struct _FmsFieldDescriptor_FixedOrder *fo = fd->descriptor.fixed_order;
  if (field_type) { *field_type = fo->field_type; }
  if (basis_type) { *basis_type = fo->basis_type; }
  if (order) { *order = fo->order; }
  return 0;
}

int FmsFieldDescriptorGetNumDofs(FmsFieldDescriptor fd, FmsInt *num_dofs) {
  if (!fd) { E_RETURN(1); }
  if (!num_dofs) { E_RETURN(2); }
  *num_dofs = fd->num_dofs;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsField functions: query interface */
/* -------------------------------------------------------------------------- */

int FmsFieldGetName(FmsField field, const char **field_name) {
  if (!field) { E_RETURN(1); }
  if (!field_name) { E_RETURN(2); }
  *field_name = field->name;
  return 0;
}

int FmsFieldGet(FmsField field, FmsFieldDescriptor *fd, FmsInt *num_vec_comp,
                FmsLayoutType *layout_type, FmsScalarType *data_type,
                const void **data) {
  if (!field) { E_RETURN(1); }
  if (fd) { *fd = field->fd; }
  if (num_vec_comp) { *num_vec_comp = field->num_vec_comp; }
  if (layout_type) { *layout_type = field->layout; }
  if (data_type) { *data_type = field->scalar_type; }
  if (data) { *data = field->data; }
  return 0;
}

int FmsFieldGetMetaData(FmsField field, FmsMetaData *mdata) {
  if (!field) { E_RETURN(1); }
  if (!mdata) { E_RETURN(2); }
  *mdata = field->mdata;
  return 0;
}


/* -------------------------------------------------------------------------- */
/* FmsMetaData functions: query interface */
/* -------------------------------------------------------------------------- */

int FmsMetaDataGetType(FmsMetaData mdata, FmsMetaDataType *type) {
  if (!mdata) { E_RETURN(1); }
  if (!type) { E_RETURN(2); }
  *type = mdata->md_type;
  return 0;
}

int FmsMetaDataGetIntegers(FmsMetaData mdata, const char **mdata_name,
                           FmsIntType *int_type, FmsInt *size,
                           const void **data) {
  if (!mdata) { E_RETURN(1); }
  if (mdata_name) { *mdata_name = mdata->name; }
  if (int_type) { *int_type = mdata->sub_type.int_type; }
  if (size) { *size = mdata->num_entries; }
  if (data) { *data = mdata->data; }
  return 0;
}

int FmsMetaDataGetScalars(FmsMetaData mdata, const char **mdata_name,
                          FmsScalarType *scal_type, FmsInt *size,
                          const void **data) {
  if (!mdata) { E_RETURN(1); }
  if (mdata_name) { *mdata_name = mdata->name; }
  if (scal_type) { *scal_type = mdata->sub_type.scalar_type; }
  if (size) { *size = mdata->num_entries; }
  if (data) { *data = mdata->data; }
  return 0;
}

int FmsMetaDataGetString(FmsMetaData mdata, const char **mdata_name,
                         const char **c_string) {
  if (!mdata) { E_RETURN(1); }
  if (mdata_name) { *mdata_name = mdata->name; }
  if (c_string) { *c_string = mdata->data; }
  return 0;
}

int FmsMetaDataGetMetaData(FmsMetaData mdata, const char **mdata_name,
                           FmsInt *size, FmsMetaData **data) {
  if (!mdata) { E_RETURN(1); }
  if (mdata_name) { *mdata_name = mdata->name; }
  if (size) { *size = mdata->num_entries; }
  if (data) { *data = mdata->data; }
  return 0;
}
