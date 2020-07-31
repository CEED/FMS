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

#ifndef FMS_HEADER
#define FMS_HEADER

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


/// Interface version constant of the form: ((major*100 + minor)*100 + patch).
enum { FMS_INTERFACE_VERSION = 100 /* v0.1 */ };

/// Type used by fms for representing and storing sizes and indices.
typedef uint64_t FmsInt;

/// TODO: dox
typedef enum {
  // signed integer types
  FMS_INT8,
  FMS_INT16,
  FMS_INT32,
  FMS_INT64,
  // unsigned integer types
  FMS_UINT8,
  FMS_UINT16,
  FMS_UINT32,
  FMS_UINT64,
  // the next constant defines the number of integer types
  // also, used as an "invalid" type
  FMS_NUM_INT_TYPES,
  /// The type of FmsInt identified as FmsIntType.
  FMS_INT_TYPE = FMS_UINT64,
  /// The type of FmsOrientation as FmsIntType.
  FMS_ORIENTATION_INT_TYPE = FMS_UINT8
} FmsIntType;

extern const size_t FmsIntTypeSize[FMS_NUM_INT_TYPES];

extern const char * const FmsIntTypeNames[FMS_NUM_INT_TYPES];

/// Get the enum representation of an int type from the string name.
int FmsGetIntTypeFromName(const char * const name, FmsIntType *type);

/// TODO: dox
/** A mesh consists of:
    * mesh domains
    * mesh components, described in terms of the domains
    * mesh tags (attributes), defined on components.

    The mesh is also assigned a partition id, e.g. an MPI rank. */
typedef struct FmsMesh_private *FmsMesh;

/// TODO: dox
/** Domains describe sets of interconnected mesh entities:
    * 0d-entities / vertices
    * 1d-entities / edges
    * 2d-entities / faces: triangles, quads
    * 3d-entities / volumes / regions: tets, hexes, wedges, pyramids.

    All entities in a domain are enumerated within the domain.

    Connections to other domains are described using shared entities. (TODO)

    Other relations, such as refinement (parent-child) and non-conforming
    (master-slave) relations can also be represented. (TODO)

    Domains are assigned a string name and an integer id. Thus, a domain is
    identified uniquely by the tuple: (name, id, partition-id) where the
    partition id is the one assigned to the containing mesh.

    The arrays storing the entity definitions use user-specified integer type.

    The domain names are stored in the mesh structure and the domains store just
    pointer to the name. */
typedef struct FmsDomain_private *FmsDomain;

/// TODO: dox
/** Components are sets of parts. A part is a set of entities of the same type,
    described as a subset of the entities of a domain. All entities in the
    component must have (i) the same dimension and (ii) specified orientation
    (relative to the entity as described in its domain).

    FIXME: parts can hold more than one type of main entities!

    In order to facilitate the definition of fields on the component, the
    following additional data can be stored in every part of the component:
    for all lower dimensional entities that are boundary to the main entities of
    the part, define an array that maps the local (to the part) indices to the
    domain-indices. These arrays plus the main array (the one describing the
    highest dimensional entities of the component) define local numberings of
    all entities inside each part. These numberings will be used to define the
    ordering of the degrees of freedom of a field. When the main entity array is
    NULL (indicating that all entities in the domain are used) then the lower
    dimensional entities will also be NULL because there is no need to have
    local numbering of the entities - the original numbering used by the domain
    can be reused.

    A component has a coordinates field (FmsField) which may be NULL.

    In addition to the parts, a component also stores relations to other
    components. A relation to a component of lower or higher dimension
    indicates a boundary or interface relation, i.e. the lower dimensional
    component describes a subset of the boundary entities of the higher
    dimensional component. A relation to another component of the same dimension
    is acceptable but has no specified meaning.

    TODO: introduce additional topological constraints in a component, e.g. to
    enforce periodicity on a non-periodic component.

    The arrays storing the entity indices into the domains use user-specified
    integer type.

    The arrays storing the entity orientations are of type FmsOrientation.

    A domain associated with a part is stored as a domain names plus a domain id
    of type FmsInt.

    The relations to other components are stored as an array of component ids of
    type FmsInt.

    Mesh tags are defined on the set of all entities in a mesh component.

    Discrete fields are defined on mesh components. Unlike tags, discrete fields
    generally associate data not only with the entities described by the mesh
    component but also with the lower-dimensional boundary entities of the
    component entities.
 */
typedef struct FmsComponent_private *FmsComponent;

/// TODO: dox
/** A tag is an array of integers describing the entities in a given component.
    Optionally, the integer tags can be assigned string descriptions.

    The array with the integer tags is ordered part-by-part within the mesh
    component.

    The associated component is stored as a component id of type FmsInt.

    A tag naturally defines a non-overlapping decomposition of the associated
    component.

    Tags can be used to store the orders (polynomial degrees) associated with
    the component entities in variable order discrete spaces.
 */
typedef struct FmsTag_private *FmsTag;

/// Constant used for initializing dimension variables.
enum { FMS_INVALID_DIM = 127 };

/// TODO: dox
/** An entity is described by a tuple of its side entity indices. For an entity
    of dimension d, its side entities are its boundary entities of dimension
    d-1.

    For FMS_TRIANGLE and FMS_QUADRILATERAL, the edges (sides), "abc"/"abcd" and
    the vertices "012"/"0123" are ordered counterclockwise, as illustrated in
    the following diagram:

        3--c--2       2
        |     |      / \
        d     b     c   b
        |     |    /     \
        0--a--1   0---a---1

    For 3D entities, the ordering of the edges inside the faces should follow
    the counterclockwise ordering when looking the the face from outside. This
    rule is followed by the choices given below.

    For FMS_TETRAHEDRON, the faces (sides), "ABCD", the edges, "abcdef", and the
    vertices, "0123" are ordered as follows:

           z=0         y=0         x=0       x+y+z=1
            2           3           3           3
           / \         / \         / \         / \
          b   c       d   e       f   d       e   f
         /  A  \     /  B  \     /  C  \     /  D  \
        1---a---0   0---a---1   2---c---0   1---b---2

    For example, vertex "0" has coordinates (x,y,z)=(0,0,0), vertex "1" has
    coordinates (1,0,0), etc.

    For FMS_HEXAHEDRON, the faces (sides), "ABCDEF", the edges, "abcdefghijkl"
    and the vertices, "01234567", are ordered as follows:

              7--g--6
             /|    /|
            / l   / k   z=0      z=1      y=0      y=1      x=0      x=1
           h  |  f  |   bottom   top      front    back     left     right
          /   3-/c--2   2--c--3  7--g--6  4--e--5  6--g--7  7--h--4  5--f--6
         /   / /   /    |     |  |     |  |     |  |     |  |     |  |     |
        4--e--5   /     b  A  d  h  B  f  i  C  j  k  D  l  l  E  i  j  F  k
        |  d  |  b      |     |  |     |  |     |  |     |  |     |  |     |
        i /   j /       1--a--0  4--e--5  0--a--1  2--c--3  3--d--0  1--b--2
        |/    |/
        0--a--1

    For example, vertex "0" has coordinates (x,y,z)=(0,0,0), vertex "6" has
    coordinates (1,1,1), etc.

    TODO: describe the side entity orderings for all ramaining entity types used
    by fms.
 */
typedef enum {
  FMS_VERTEX,
  FMS_EDGE,
  FMS_TRIANGLE,
  FMS_QUADRILATERAL,
  FMS_TETRAHEDRON,
  FMS_HEXAHEDRON,
  FMS_WEDGE, // triangular prism (deformed)
  FMS_PYRAMID,
  // the next constant defines the number of entity types
  FMS_NUM_ENTITY_TYPES
} FmsEntityType;

/// String representations of each entity type.
extern const char * const FmsEntityTypeNames[FMS_NUM_ENTITY_TYPES];

int FmsGetEntityTypeFromName(const char * const name, FmsEntityType *ent_type);

/// Dimensions of the entity types.
extern const FmsInt FmsEntityDim[FMS_NUM_ENTITY_TYPES];

/// Number of sides of the entity types.
/** The sides of an entity of dimension d are its boundary entities of dimension
    d-1. */
extern const FmsInt FmsEntityNumSides[FMS_NUM_ENTITY_TYPES];

/// Number of vertices of the entity types.
extern const FmsInt FmsEntityNumVerts[FMS_NUM_ENTITY_TYPES];

/** @brief A type describing user-defined local orderings of the side entities
    defining an enitity. */
/** User-ordering to fms-ordering is performed as follows: let u=[u_0,...,u_n]
    be the user-ordered indices for one of the entity types, say t, which has
    reordering array ra (of size n+1), i.e. ra is the pointer given by entry t
    from a variable of type FmsEntityReordering (which is an array of pointers).
    Then the fms-ordered indices, v=[v_0,...,v_n], are computed as v[i] =
    u[ra[i]], i=0,...,n. Conversely, if v is given, then u is computed as
    u[ra[i]] = v[i], i=0,...,n. In other words, ra[fms-index] = user-index. If a
    pointer ra is NULL, then no reordering is performed.

    This type can be used to describe orderings of other types of boundary
    entities, e.g. reorderings for the vertices.
 */
typedef const int *FmsEntityReordering[FMS_NUM_ENTITY_TYPES];

/// A type used by fms to represent and store entity orientations.
/** Orientations identify specific permutations to the d-1 dimensional boundary
    entities that describe a d-dimensional entity. Orientations are always
    relative to a specific entity described as a tuple of d-1 dimensional
    entities. For example, if the triangle A, given as A=(a,b,c) where a, b, and
    c are edges (i.e. edge indices) then we say that the triangle B=(b,c,a) has
    orientation 2 relative to A, see the definitions below.

    In general, the orientation of an entity B relative to an entity A is
    different from the orientation of A relative to B. These two orientations
    are inversees of each other: mathematically the set of all orientation of an
    entity type form a group.

    Orientation 0 always represents the identity permutation.

    For FMS_EDGE, the orientations of an edge "ab" are defined as:

        0: ab   1: ba

    For FMS_TRIANGLE, the orientations of a triangle "abc" are defined as:

        0: abc   1: acb   2: bca   3: bac   4: cab   5: cba

    The inverses of the orientations `012345` are the orientations `014325`,
    respectively, i.e. the orientations `2` and `4` are inverses of each other
    while all other orientations are inverses of themselves.

    For FMS_QUADRILATERAL, the orientations of a quad "abcd" are defined as:

        0: abcd   1: adcb   2: bcda   3: badb
        4: cdab   5: cbad   6: dabc   7: dcba

    The inverses of the orientations `01234567` are the orientations `01634527`,
    respectively, i.e. `2` and `6` are inverses of each other while all other
    orientations are inverses of themselves.
 */
typedef uint8_t FmsOrientation;

enum { FMS_ORIENTATION_UNKNOWN = 255 };

/// TODO: dox
typedef enum {
  FMS_FLOAT,
  FMS_DOUBLE,
  FMS_COMPLEX_FLOAT,
  FMS_COMPLEX_DOUBLE,
  // next constant gives the number of scalar types
  FMS_NUM_SCALAR_TYPES
} FmsScalarType;

extern const size_t FmsScalarTypeSize[FMS_NUM_SCALAR_TYPES];

extern const char * const FmsScalarTypeNames[FMS_NUM_SCALAR_TYPES];

/// Get the enum representation of an int type from the string name.
int FmsGetScalarTypeFromName(const char * const name, FmsScalarType *type);

/// TODO: dox
typedef enum {
  FMS_FIXED_ORDER
} FmsFieldDescriptorType;

/// TODO: dox
typedef enum {
  FMS_CONTINUOUS,    // H1-conforming
  FMS_DISCONTINUOUS, // "L2-conforming"
  FMS_DISCONTINUOUS_WEIGHTED, // "L2-conforming" volume weighted
  FMS_HCURL,         // H(curl)-conforming
  FMS_HDIV           // H(div)-conforming
} FmsFieldType;

/// TODO: dox
typedef enum {
  FMS_NODAL_GAUSS_OPEN,
  FMS_NODAL_GAUSS_CLOSED,
  FMS_NODAL_CHEBYSHEV_OPEN,
  FMS_NODAL_CHEBYSHEV_CLOSED,
  FMS_NODAL_UNIFORM_OPEN,
  FMS_NODAL_UNIFORM_CLOSED,
  FMS_POSITIVE
} FmsBasisType;

/// TODO: dox
typedef enum {
  FMS_BY_NODES, // "XX..., YY..., ZZ..."
  FMS_BY_VDIM   // "XYZ, XYZ, ..."
} FmsLayoutType;

/// TODO: dox
typedef enum {
  FMS_INTEGER,
  FMS_SCALAR,
  FMS_STRING,
  FMS_META_DATA,
  FMS_NUM_METADATA_TYPES
} FmsMetaDataType;

extern const char * const FmsMetaDataTypeNames[FMS_NUM_METADATA_TYPES];

int FmsGetMetaDataTypeFromName(const char * const name, FmsMetaDataType *type);

/// TODO: dox
/** A meta-data structure contains:
    * a meta-data type (FmsMetaDataType)
    * a meta-data subtype, e.g. FmsIntType, FmsScalarType; TODO: for FMS_STRING
      type define FmsEncoding subtype
    * a meta-data name (const char *)
    * number of entries in the data array (FmsInt); the type of the entries in
      the data array is based on the meta-data type and subtype, when subtype is
      applicable; in the case of FMS_STRING, the entries are of type char (or
      some other type depending on the FmsEncoding subtype when introduced)
    * a data array (void *)
 */
typedef struct FmsMetaData_private *FmsMetaData;

/// TODO: dox
/** A field-descriptor structure contains:
    * a field-descriptor name
    * an associated mesh component (FmsComponent)
    * a descriptor-type (FmsFieldDescriptorType)
    * if descriptor-type == "fixed-order":
      - a field type (FmsFieldType)
      - a basis-type (FmsBasisType)
      - an order
    * total number of DOFs

    The degrees of freedom are ordered part-by-part of the associated mesh
    component. Within each part, the dofs are ordered as follows:
    * first, all DOFs on all vertices ordered according to the local (to the
      part) enumeration of the vertices
    * second, all DOFs on all edges ordered according to the local (to the part)
      enumeration of the edges
    * similarly, continue adding DOFs for all remaining entity types in the
      order used by FmsEntityType.
 */
typedef struct FmsFieldDescriptor_private *FmsFieldDescriptor;

/// TODO: dox
/** A field structure contains:
    * a field name
    * meta-data - e.g. time (FmsMetaData, can be NULL)
    * a field-descriptor (FmsFieldDescriptor)
    * number of vector components (FmsInt)
    * a layout type, i.e. ordering for the vector components (FmsLayoutType)
    * a scalar type (FmsScalarType)
    * a data array of the given scalar type.
 */
typedef struct FmsField_private *FmsField;

/// TODO: dox
/** A data collection structure contains:
   * a name (const char *)
   * meta-data (FmsMetaData, can be NULL)
   * a mesh (FmsMesh)
   * number of field-descriptors (FmsInt)
   * an array of field-descriptors (FmsFieldDescriptor *)
   * number of fields (FmsInt)
   * an array of fields (FmsField *)
 */
typedef struct FmsDataCollection_private *FmsDataCollection;


//
//   General functions
//

/// Get the interface version.
/** Version is of the form: ((major*100 + minor)*100 + patch). */
int FmsGetInterfaceVersion(FmsInt *version);


//
//   Construction interface
//

/// Allocate a mesh structure and initialize it to be empty.
/** After this call, the content of the mesh should be described by calling the
    functions FmsMeshAdd* and, if the mesh is partitioned,
    FmsMeshSetPartitionId.

    To complete the construction of the mesh structure, the function
    FmsMeshFinalize must be called. */
int FmsMeshConstruct(FmsMesh *mesh);

/// TODO: dox
int FmsMeshFinalize(FmsMesh mesh);

/// TODO: dox
int FmsMeshValidate(FmsMesh mesh);

/// TODO: dox
int FmsMeshDestroy(FmsMesh *mesh);

/// TODO: dox
int FmsMeshSetPartitionId(FmsMesh mesh, FmsInt partition_id,
                          FmsInt num_partitions);

/** @brief Allocates an array of domains sharing the same name, initializes
    them, and returns a pointer to the array. */
/** After this call, the domains should be described via the FmsDomainSet* and
    FmsDomainAdd* functions. */
int FmsMeshAddDomains(FmsMesh mesh, const char *domain_name, FmsInt num_domains,
                      FmsDomain **domains);

/** @brief Add a new empty component to the mesh with the given name and return
    the new component in @a comp. */
/** After this call, the component should be described via the FmsComponentAdd*
    functions. */
int FmsMeshAddComponent(FmsMesh mesh, const char *comp_name,
                        FmsComponent *comp);

/// TODO: dox
/** @brief Add a new tag to the mesh with the given name and return the new tag
    in @a tag. */
/** The tag should be set via the FmsTagSet* functions. */
int FmsMeshAddTag(FmsMesh mesh, const char *tag_name, FmsTag *tag);

/// Set the number of vertices in a domain.
int FmsDomainSetNumVertices(FmsDomain domain, FmsInt num_verts);

/** @brief Allocates memory for the specified entities. Destroys any existing
    entities of the given type. */
/** The value of @a num_ents can be an upper bound of the actual number of
    entities. */
int FmsDomainSetNumEntities(FmsDomain domain, FmsEntityType type,
                            FmsIntType id_store_type, FmsInt num_ents);

/// TODO: dox
/** This function allows the internal entities array to be set directly by the
    user. */
int FmsDomainGetEntitiesArray(FmsDomain domain, FmsEntityType type,
                              void **ents);

/// TODO: dox
/** This function allows the internal entity side orientations array to be set
    directly by the user. */
int FmsDomainGetOrientationsArray(FmsDomain domain, FmsEntityType type,
                                  FmsOrientation **side_orients);

/// TODO: dox
/** Entities of dimension d are described by tuples of its side entities (entity
    indices), or in other words, by a list of its boundary entities of dimension
    d-1.

    The content of the array @a ents is copied internally after applying the
    side entity reordering. */
int FmsDomainAddEntities(FmsDomain domain, FmsEntityType type,
                         FmsEntityReordering reordering, FmsIntType ent_id_type,
                         const void *ents, FmsInt num_ents);

/// TODO: dox
int FmsDomainAddOrientation(FmsDomain domain, FmsEntityType type, FmsInt ent_id,
                            FmsInt loc_side_id, FmsOrientation side_orient);

/// TODO: dox
/** All entities of the highest dimension in the domain are added to the
    component. */
int FmsComponentAddDomain(FmsComponent comp, FmsDomain domain);

/// TODO: dox
int FmsComponentAddPart(FmsComponent comp, FmsDomain domain, FmsInt *part_id);

/// TODO: dox
/** The pointer @a ents can be NULL to signify that all entities in the domain
    of the given type are being added.

    The pointer @a ent_orients can be NULL to signify that all orientations are
    0, i.e. the identity permutation.

    The arrays @a ents and @a ent_orients (if not NULL) have the same size, @a
    num_ents.

    The content of the arrays @a ents and @a ent_orients is copied internally
    after applying the @a inv_orient_map to the second array.

    TODO: inv_orient_map[user-orientation] -> fms-orientation. */
int FmsComponentAddPartEntities(FmsComponent comp, FmsInt part_id,
                                FmsEntityType type, FmsIntType id_store_type,
                                FmsIntType ent_id_type, FmsIntType orient_type,
                                const FmsOrientation *inv_orient_map,
                                const void *ents, const void *ent_orients,
                                FmsInt num_ents);

/// TODO: dox
/// Similar to FmsComponentAddPartEntities() but for lower dimensional entities.
int FmsComponentAddPartSubEntities(FmsComponent comp, FmsInt part_id,
                                   FmsEntityType type, FmsIntType id_store_type,
                                   FmsIntType ent_id_type, const void *ents,
                                   FmsInt num_ents);

/// Describe a relation from one component to another.
int FmsComponentAddRelation(FmsComponent comp, FmsInt other_comp_id);

/// Set the coordinates field of a component.
int FmsComponentSetCoordinates(FmsComponent comp, FmsField coords);

/// TODO: dox
int FmsTagSetComponent(FmsTag tag, FmsComponent comp);

/// TODO: dox
/** The array @a ent_tags is converted and copied internally. The value of
    @a num_ents must be the same as the number of main entities in the
    associated mesh component. */
int FmsTagSet(FmsTag tag, FmsIntType stored_tag_type, FmsIntType input_tag_type,
              const void *ent_tags, FmsInt num_ents);

/// TODO: dox
/** The pointer @a num_ents can be NULL. */
int FmsTagAllocate(FmsTag tag, FmsIntType stored_tag_type, void **ent_tags,
                   FmsInt *num_ents);

/// TODO: dox
/** The arrays @a tags and @a tag_descr have the same size, @a num_tags. The
    tags and their descriptions are copied internally. The @a tag_type is the
    type of the input array @a tags. */
int FmsTagAddDescriptions(FmsTag tag, FmsIntType tag_type, const void *tags,
                          const char *const *tag_descr, FmsInt num_tags);

/// TODO: dox
/// The new object, @a dc, assumes ownership of the @a mesh.
int FmsDataCollectionCreate(FmsMesh mesh, const char *dc_name,
                            FmsDataCollection *dc);

/// Destroy a data collection object.
/** Destroys the mesh, all field-descriptors, and all fields stored by the data
    collection. */
int FmsDataCollectionDestroy(FmsDataCollection *dc);

/// TODO: dox
int FmsDataCollectionAddFieldDescriptor(FmsDataCollection dc,
                                        const char *fd_name,
                                        FmsFieldDescriptor *fd);

/// TODO: dox
int FmsDataCollectionAddField(FmsDataCollection dc, const char *field_name,
                              FmsField *field);

/** @brief Make sure the meta-data structure associated with a data collection
    is allocated and return it in @a mdata. */
int FmsDataCollectionAttachMetaData(FmsDataCollection dc, FmsMetaData *mdata);

/// TODO: dox
int FmsFieldDescriptorSetComponent(FmsFieldDescriptor fd, FmsComponent comp);

/// TODO: dox
int FmsFieldDescriptorSetFixedOrder(FmsFieldDescriptor fd,
                                    FmsFieldType field_type,
                                    FmsBasisType basis_type, FmsInt order);

/// TODO: dox
/// The size of the @a data array is @a num_vec_comp times the number of DOFs
/// as defined by the field descriptor @a fd.
int FmsFieldSet(FmsField field, FmsFieldDescriptor fd, FmsInt num_vec_comp,
                FmsLayoutType layout_type, FmsScalarType data_type,
                const void *data);

/** @brief Make sure the meta-data structure associated with a field is
    allocated and return it in @a mdata. */
int FmsFieldAttachMetaData(FmsField field, FmsMetaData *mdata);

/// Set the contents of a meta-data structure to store an array of integers.
/** Returns the uninitialized meta-data array in @a data. */
int FmsMetaDataSetIntegers(FmsMetaData mdata, const char *mdata_name,
                           FmsIntType int_type, FmsInt size, void **data);

/// Set the contents of a meta-data structure to store an array of scalars.
/** Returns the uninitialized meta-data array in @a data. */
int FmsMetaDataSetScalars(FmsMetaData mdata, const char *mdata_name,
                          FmsScalarType scal_type, FmsInt size, void **data);

/// Set the contents of a meta-data structure to store a c-string.
int FmsMetaDataSetString(FmsMetaData mdata, const char *mdata_name,
                         const char *c_string);

/** @brief Set the contents of a meta-data structure to store an array of
    meta-data structures. */
/** Returns the meta-data array in @a data. Each meta-data entry in the array is
    initialized to be empty. */
int FmsMetaDataSetMetaData(FmsMetaData mdata, const char *mdata_name,
                           FmsInt size, FmsMetaData **data);

/// Clear the contents of @a mdata without destroying the object.
int FmsMetaDataClear(FmsMetaData mdata);

/// Destroy the object that @a *mdata refers to and set @a *mdata to NULL.
int FmsMetaDataDestroy(FmsMetaData *mdata);


//
//   Query interface
//

/// TODO: dox
int FmsMeshGetPartitionId(FmsMesh mesh, FmsInt *partition_id,
                          FmsInt *num_partitions);

/// TODO: dox
int FmsMeshGetNumDomainNames(FmsMesh mesh, FmsInt *num_domain_names);

/// TODO: dox
int FmsMeshGetDomains(FmsMesh mesh, FmsInt domain_name_id,
                      const char **domain_name, FmsInt *num_domains,
                      FmsDomain **domains);

/// TODO: dox
int FmsMeshGetDomainsByName(FmsMesh mesh, const char *domain_name,
                            FmsInt *num_domains, FmsDomain **domains);

/// TODO: dox
int FmsMeshGetNumComponents(FmsMesh mesh, FmsInt *num_comp);

/// TODO: dox
int FmsMeshGetComponent(FmsMesh mesh, FmsInt comp_id, FmsComponent *comp);

/// TODO: dox
int FmsMeshGetNumTags(FmsMesh mesh, FmsInt *num_tags);

/// TODO: dox
int FmsMeshGetTag(FmsMesh mesh, FmsInt tag_id, FmsTag *tag);

/// Return the name and id of a mesh domain.
int FmsDomainGetName(FmsDomain domain, const char **domain_name,
                     FmsInt *domain_id);

/// Return the highest dimension of an entry in a mesh domain.
int FmsDomainGetDimension(FmsDomain domain, FmsInt *dim);

/// Get the number of vertices in a domain.
int FmsDomainGetNumVertices(FmsDomain domain, FmsInt *num_verts);

/// Get the number of entities of a given type in a domain.
int FmsDomainGetNumEntities(FmsDomain domain, FmsEntityType type,
                            FmsInt *num_ents);

/// TODO: dox
/// No copy, read only access to all entity definitions.
int FmsDomainGetAllEntities(FmsDomain domain, FmsEntityType type,
                            FmsIntType *ent_id_type, const void **ents,
                            FmsInt *num_ents);

/// TODO: dox
/// Extract a subset of the entity definitions.
int FmsDomainGetEntities(FmsDomain domain, FmsEntityType type,
                         FmsEntityReordering reordering, FmsIntType ent_id_type,
                         FmsInt first_ent, void *ents, FmsInt num_ents);

/// TODO: dox
/// Extract a subset of the entity definitions in terms of their vertices.
int FmsDomainGetEntitiesVerts(FmsDomain domain, FmsEntityType type,
                              FmsEntityReordering vert_reordering,
                              FmsIntType ent_id_type, FmsInt first_ent,
                              void *ents_verts, FmsInt num_ents);

/// TODO: dox
/// No copy, read only access to all entity side orientations.
int FmsDomainGetAllOrientations(FmsDomain domain, FmsEntityType type,
                                const FmsOrientation **side_orients,
                                FmsInt *num_ents);

/// Return the name of a mesh component.
int FmsComponentGetName(FmsComponent comp, const char **comp_name);

/// Return the dimension of a mesh component.
int FmsComponentGetDimension(FmsComponent comp, FmsInt *dim);

/// TODO: dox
int FmsComponentGetNumParts(FmsComponent comp, FmsInt *num_parts);

/// TODO: dox
/** No copy, read only access to the description of part @a part_id in the
    component.

    If the pointer returned in @a ents is NULL, then all entities in the domain
    of the specified type are used.

    If the pointer returned in @a ent_orients is NULL, then all orientations are
    0, i.e. the identity permutation. */
int FmsComponentGetPart(FmsComponent comp, FmsInt part_id, FmsEntityType type,
                        FmsDomain *domain, FmsIntType *ent_id_type,
                        const void **ents, const FmsOrientation **ent_orients,
                        FmsInt *num_ents);

/// TODO: dox
/// Similar to FmsComponentGetPart() but for lower dimensional entities.
int FmsComponentGetPartSubEntities(FmsComponent comp, FmsInt part_id,
                                   FmsEntityType type, FmsIntType *ent_id_type,
                                   const void **ents, FmsInt *num_ents);

/// Return the total number of main entities across all parts in the component.
int FmsComponentGetNumEntities(FmsComponent comp, FmsInt *num_ents);

/// TODO: dox
int FmsComponentGetRelations(FmsComponent comp, const FmsInt **rel_comps,
                             FmsInt *num_rel_comps);

/// Return the coordinates field of a component.
/** Note that the returned FmsField, @a *coords, can be NULL. */
int FmsComponentGetCoordinates(FmsComponent comp, FmsField *coords);

/// TODO: dox
int FmsTagGetName(FmsTag tag, const char **tag_name);

/// TODO: dox
int FmsTagGetComponent(FmsTag tag, FmsComponent *comp);

/// TODO: dox
/// No copy, read only access to the entity-tag array.
int FmsTagGet(FmsTag tag, FmsIntType *tag_type, const void **ent_tags,
              FmsInt *num_ents);

/// TODO: dox
/** The arrays @a tags and @a tag_descr have the same size, @a num_tags.

    The @a tag_type is the type of the entries in the @a tags array. */
int FmsTagGetDescriptions(FmsTag tag, FmsIntType *tag_type, const void **tags,
                          const char *const **tag_descr, FmsInt *num_tags);

/// TODO: dox
/**
@brief Get the name of the data collection.
@param dc the data collection.
@param[out] name A char pointer that will be set to point to the name string.
           The name is owned by the data collection so it should not be 
           freed or modified.
@return 0 on success, non-zero otherwise.
*/
int FmsDataCollectionGetName(FmsDataCollection dc, char **name);

/// TODO: dox
int FmsDataCollectionGetMesh(FmsDataCollection dc, FmsMesh *mesh);

/// TODO: dox
int FmsDataCollectionGetFieldDescriptors(FmsDataCollection dc,
    FmsFieldDescriptor **fds, FmsInt *num_fds);

/// TODO: dox
int FmsDataCollectionGetFields(FmsDataCollection dc, FmsField **fields,
                               FmsInt *num_fields);

/// TODO: dox
int FmsDataCollectionGetMetaData(FmsDataCollection dc, FmsMetaData *mdata);

/// TODO: dox
int FmsFieldDescriptorGetName(FmsFieldDescriptor fd, const char **fd_name);

/// TODO: dox
int FmsFieldDescriptorGetComponent(FmsFieldDescriptor fd, FmsComponent *comp);

/// TODO: dox
int FmsFieldDescriptorGetType(FmsFieldDescriptor fd,
                              FmsFieldDescriptorType *fd_type);

/// TODO: dox
int FmsFieldDescriptorGetFixedOrder(FmsFieldDescriptor fd,
                                    FmsFieldType *field_type,
                                    FmsBasisType *basis_type, FmsInt *order);

/// TODO: dox
int FmsFieldDescriptorGetNumDofs(FmsFieldDescriptor fd, FmsInt *num_dofs);

/// TODO: dox
int FmsFieldGetName(FmsField field, const char **field_name);

/// TODO: dox
int FmsFieldGet(FmsField field, FmsFieldDescriptor *fd, FmsInt *num_vec_comp,
                FmsLayoutType *layout_type, FmsScalarType *data_type,
                const void **data);

/// TODO: dox
int FmsFieldGetMetaData(FmsField field, FmsMetaData *mdata);

/// TODO: dox
/** If the @a mdata structure is empty, returns an error code of 1. */
int FmsMetaDataGetType(FmsMetaData mdata, FmsMetaDataType *type);

/// Get the contents of a meta-data structure that stores an array of integers.
/** Returns a pointer to the meta-data array in @a data. */
int FmsMetaDataGetIntegers(FmsMetaData mdata, const char **mdata_name,
                           FmsIntType *int_type, FmsInt *size,
                           const void **data);

/// Get the contents of a meta-data structure that stores an array of scalars.
/** Returns a pointer to the meta-data array in @a data. */
int FmsMetaDataGetScalars(FmsMetaData mdata, const char **mdata_name,
                          FmsScalarType *scal_type, FmsInt *size,
                          const void **data);

/// Get the contents of a meta-data structure that stores a c-string.
int FmsMetaDataGetString(FmsMetaData mdata, const char **mdata_name,
                         const char **c_string);

/** @brief Get the contents of a meta-data structure that stores an array of
    meta-data structures. */
/** Returns a pointer to the meta-data array in @a data. */
int FmsMetaDataGetMetaData(FmsMetaData mdata, const char **mdata_name,
                           FmsInt *size, FmsMetaData **data);


#ifdef __cplusplus
} // extern "C"
#endif

#endif // FMS_HEADER
