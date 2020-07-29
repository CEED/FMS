/*
 Copyright (c) 2020, Lawrence Livermore National Security, LLC. Produced at
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

#include <fmsio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#ifdef FMS_HAVE_CONDUIT
#include <conduit.h>
#include <conduit_relay.h>
#endif

#define FREE(PTR) if((PTR) != NULL) free(PTR)

// Feature: If you set this to 1 then the last line is going to have 2
#define FMS_ELE_PER_LINE 3u

#if UINT_MAX == ULONG_MAX
#define FMS_LD "%lld"
#define FMS_LU "%llu"
#else
#define FMS_LD "%ld"
#define FMS_LU "%lu"
#endif

#define FOR_EACH_INT_TYPE(macro)            \
    macro(FMS_INT8, int8_t,     "%d")       \
    macro(FMS_INT16, int16_t,   "%d")       \
    macro(FMS_INT32, int32_t,   "%d")       \
    macro(FMS_INT64, int64_t,   FMS_LD)    \
    macro(FMS_UINT8, uint8_t,   "%u")       \
    macro(FMS_UINT16, uint16_t, "%u")       \
    macro(FMS_UINT32, uint32_t, "%u")       \
    macro(FMS_UINT64, uint64_t, FMS_LU)

#define FOR_EACH_SCALAR_TYPE(macro)         \
    macro(FMS_FLOAT, float,   "%f")         \
    macro(FMS_DOUBLE, double, "%f")         \
    macro(FMS_COMPLEX_FLOAT, float, "%f")   \
    macro(FMS_COMPLEX_DOUBLE, double, "%f")

/**
@note borrowed from fms.c for a bit.
*/
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

/**
Data structures.
*/
typedef struct
{
#ifdef FMS_HAVE_CONDUIT
    /* Conduit C-API*/
    conduit_node *root;

#ifdef MEANDERING_THOUGHTS
    conduit_node *stack[100];
    int           stack_top;
#endif
    char *filename;
    char *protocol;
    int   writing;
    
#endif
    FILE *fp;
} FmsIOContext;

typedef struct
{
    /* Function pointers */
    int (*open)(FmsIOContext *ctx, const char *filename, const char *mode);
    int (*close)(FmsIOContext *ctx);

#ifdef MEANDERING_THOUGHTS
    /* Q: would it be helpful to have some kind of begin/end object and list functions?

      I implemented some of this and then had some second thoughts.
      Leaving it here for a little while.
     */
    int (*begin_list)(FmsIOContext *ctx, const char *path, size_t n);
    int (*end_list)(FmsIOContext *ctx);
    int (*begin_list_item)(FmsIOContext *ctx);
    int (*end_list_item)(FmsIOContext *ctx);
#endif

    int (*add_int)(FmsIOContext *ctx, const char *path, FmsInt value);
    int (*add_int_array)(FmsIOContext *ctx, const char *path, const FmsInt *values, FmsInt n);
    int (*add_typed_int_array)(FmsIOContext *ctx, const char *path, FmsIntType type, const void *values, FmsInt n);
    int (*add_float)(FmsIOContext *ctx, const char *path, float value);
    /* int (*add_float_array)(FmsIOContext *ctx, const char *path, const float *values, FmsInt n); */
    int (*add_double)(FmsIOContext *ctx, const char *path, double value);
    /* int (*add_double_array)(FmsIOContext *ctx, const char *path, const double *value, FmsInt n); */
    int (*add_scalar_array)(FmsIOContext *ctx, const char *path, FmsScalarType type, const void *data, FmsInt n);
    int (*add_string)(FmsIOContext *ctx, const char *path, const char *value);
    int (*add_string_array)(FmsIOContext *ctx, const char *path, const char **value, FmsInt n);

    /* TODO: functions for reading path+value*/
    int (*has_path)(FmsIOContext *ctx, const char *path);
    int (*get_int)(FmsIOContext *ctx, const char *path, FmsInt *value);
    int (*get_typed_int_array)(FmsIOContext *ctx, const char *path, FmsIntType *type, void **values, FmsInt *n);
    int (*get_float)(FmsIOContext *ctx, const char *path, float *value);
    int (*get_double)(FmsIOContext *ctx, const char *path, double *value);
    int (*get_scalar_array)(FmsIOContext *ctx, const char *path, FmsScalarType *type, void **values, FmsInt *n);
    int (*get_string)(FmsIOContext *ctx, const char *path, const char **value);
    int (*get_string_array)(FmsIOContext *ctx, const char *path, const char ***value, FmsInt *n);
    /*...*/

} FmsIOFunctions;

// Struct used for building metadata
typedef struct {
    FmsMetaDataType mdtype;
    union {
        FmsIntType i_type;
        FmsScalarType s_type;
    };
    const char *name;
    FmsInt size;    
    void *data;
} FmsIOMetaDataInfo;

// Struct used for building a FieldDescriptor on the read end
typedef struct {
    const char *name;
    const char *component_name;
    FmsFieldDescriptorType type;
    FmsInt fixed_order[3];
    FmsInt num_dofs;
} FmsIOFieldDescriptorInfo;

// Struct used for building a Field on the read end
typedef struct {
    const char *name;
    FmsLayoutType layout;
    FmsInt num_vec_comps;
    const char *fd_name;
    FmsInt data_size;
    FmsScalarType data_type;
    void *data; // Currently - should free
} FmsIOFieldInfo;

// Struct used to define entities in a domain
typedef struct {
    FmsEntityType ent_type;
    FmsInt nents;
    FmsInt size;
    FmsIntType type;
    void *values; // Currently - should free
} FmsIOEntityInfo;

// Struct used to build a domain
typedef struct {
    FmsInt dim; // Implies length of entities array
    FmsInt nverts;
    FmsIOEntityInfo *entities;
} FmsIODomainInfo;

// Struct used for holding domain name info
typedef struct {
    const char *name;
    FmsInt ndomains;
    FmsIODomainInfo *domains;
} FmsIODomainNameInfo;

typedef struct {
    const char *domain_name;
    FmsInt domain_id;
    FmsInt full_domain; // Boolean
    FmsEntityType part_ent_type;
    FmsIOEntityInfo entities[FMS_NUM_ENTITY_TYPES];
} FmsIOComponentPartInfo;

// Struct used for building components
typedef struct {
    const char *name;
    FmsInt dim;
    FmsInt nents;
    const char *coordinates_name;
    FmsInt nparts;
    FmsIOComponentPartInfo *parts;
    FmsInt relation_size;
    FmsIntType relation_type;
    void *relation_values;
} FmsIOComponentInfo;

// Struct used for building tags
typedef struct {
    const char *name;
    const char *comp_name;
    FmsIntType type;
    void *values;
    FmsInt size;
    FmsInt descr_size;
    FmsInt *tag_values; // For descriptions, TODO: Support unsigned
    const char **descriptions;
} FmsIOTagInfo;

// Struct used for building a Mesh
typedef struct {
    FmsInt partition_info[2];
    FmsInt ndomain_names;
    FmsInt ncomponents;
    FmsInt ntags;
    FmsIODomainNameInfo *domain_names;
    FmsIOComponentInfo *components;
    FmsIOTagInfo *tags;

} FmsIOMeshInfo;

// Struct used for building a datacollection
typedef struct {
    const char *name;
    FmsInt nfds;
    FmsIOFieldDescriptorInfo *fds;
    FmsInt nfields;
    FmsIOFieldInfo *fields;
    FmsIOMeshInfo mesh_info;
    FmsIOMetaDataInfo *md;
} FmsIODataCollectionInfo;

static char *
join_keys(const char *k1, const char *k2)
{
    size_t len = 0;
    char *str = NULL;

    len = strlen(k1) + 1/*slash*/ + strlen(k2) + 1;
    str = (char *)malloc(sizeof(char) * len);
    if(str)
        sprintf(str, "%s/%s", k1, k2);
    return str;
}

/*****************************************************************************
*** FMS I/O Functions for ASCII
*****************************************************************************/
static int
FmsIOOpen(FmsIOContext *ctx, const char *filename, const char *mode)
{
    if(!ctx) E_RETURN(1);
    if(!filename) E_RETURN(2);
    if(!mode) E_RETURN(3);
    ctx->fp = fopen(filename, mode);
    if(!ctx->fp) E_RETURN(4);
    return 0;
}

static int
FmsIOClose(FmsIOContext *ctx)
{
    if(ctx == NULL) E_RETURN(1);
    fclose(ctx->fp);
    ctx->fp = NULL;
    return 0;
}

/**
NOTE: These functions are where we can control what the ASCII file format looks like
*/
static int
FmsIOAddInt(FmsIOContext *ctx, const char *path, FmsInt value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    fprintf(ctx->fp, "%s: " FMS_LU "\n", path, value);
    return 0;
}

static int
FmsIOAddIntArray(FmsIOContext *ctx, const char *path, const FmsInt *values, FmsInt n)
{
    FmsInt i, j;
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    fprintf(ctx->fp, "%s/Size: " FMS_LU "\n", path, n);
    fprintf(ctx->fp, "%s/Type: %s\n", path, FmsIntTypeNames[FMS_UINT64]);

#if 1
    /* Should we make it YAML-like?*/
    fprintf(ctx->fp, "%s/Values: [", path);
    for(i = 0; i < n; ++i)
    {
        for(j = 0; j < FMS_ELE_PER_LINE && i < n; j++, i++)
        {
            fprintf(ctx->fp, FMS_LU, values[i]);
            if(i < n-1)
                fprintf(ctx->fp, ", ");
        }
        if(i < n-1)
            fprintf(ctx->fp, "\n");
    }
    fprintf(ctx->fp, "]\n");
#else
    /* Or should we make it a little simpler to read (write len)?*/
    fprintf(ctx->fp, "%s: %llu ", path, n);
    for(i = 0; i < n; ++i)
        fprintf(ctx->fp, "%llu ", values[i]);
    fprintf(ctx->fp, "\n");
#endif
    return 0;
}


/**
NOTE: This function is needed since for metadata we're given void* data with an enum that describes the type.
      In our ascii format, do we want to consider adding the types to make it unambiguous when we read back?

      If we added the type as a comment at the end, it would still be valid YAML.

      key1: value                             # string
      key2/path/to/something: "blah"          # string
      key2/path/to/other: 1.2345              # float
      key2/path/to/other2: [0,1,2,3,4,5,6]    # int
*/
static int
FmsIOAddTypedIntArray(FmsIOContext *ctx, const char *path, FmsIntType type, const void *values, FmsInt n)
{
    /* Shall we stick some type information into the ASCII to preserve the intended IntType? */
    size_t i, j;
    int retval = 0;
    const unsigned int epl = FMS_ELE_PER_LINE;

    fprintf(ctx->fp, "%s/Size: " FMS_LU "\n", path, n);
    fprintf(ctx->fp, "%s/Type: %s\n", path, FmsIntTypeNames[type]);

#define PRINT_MACRO(typename, T, format)    \
    do {                                    \
    const T *v = (const T *)(values);       \
    fprintf(ctx->fp, "%s/Values: [", path); \
    for(i = 0; i < n;) {                    \
        for(j = 0; j < epl && i < n; j++, i++) { \
            fprintf(ctx->fp, format, v[i]); \
            if(i < n-1)                     \
                fprintf(ctx->fp, ", ");     \
        }                                   \
        if(i < n-1)                         \
            fprintf(ctx->fp, "\n");         \
    }                                       \
    fprintf(ctx->fp, "]\n");                \
    } while(0)

    switch(type)
    {
    #define CASES(typename, T, format)          \
        case typename:                          \
            PRINT_MACRO(typename, T, format);   \
            break;
        FOR_EACH_INT_TYPE(CASES)
    #undef CASES
        default:
            E_RETURN(1);
            break;
    }
#undef PRINT_MACRO
    return retval;
}

static int
FmsIOAddFloat(FmsIOContext *ctx, const char *path, float value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    fprintf(ctx->fp, "%s: %f\n", path, value);
    return 0;
}

static int
FmsIOAddDouble(FmsIOContext *ctx, const char *path, double value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    fprintf(ctx->fp, "%s: %f\n", path, value);
    return 0;
}

static int
FmsIOAddScalarArray(FmsIOContext *ctx, const char *path, FmsScalarType type, const void *values, FmsInt n)
{
    /* Shall we stick some type information into the ASCII to preserve the intended IntType? */
    size_t i, j;
    int retval = 0;
    const unsigned int epl = FMS_ELE_PER_LINE;

    if(type == FMS_COMPLEX_FLOAT || type == FMS_COMPLEX_FLOAT)
        n = n * 2;

    fprintf(ctx->fp, "%s/Size: " FMS_LU "\n", path, n);
    fprintf(ctx->fp, "%s/Type: %s\n", path, FmsScalarTypeNames[type]);

#define PRINT_MACRO(typename, T, format)    \
    do {                                    \
    const T *v = (const T *)(values);       \
    fprintf(ctx->fp, "%s/Values: [", path); \
    for(i = 0; i < n;) {                    \
        for(j = 0; j < epl && i < n; j++, i++) { \
            fprintf(ctx->fp, format, v[i]); \
            if(i < n-1)                     \
                fprintf(ctx->fp, ", ");     \
        }                                   \
        if(i < n-1)                         \
            fprintf(ctx->fp, "\n");         \
    }                                       \
    fprintf(ctx->fp, "]\n");                \
    } while(0)

    switch(type)
    {
    #define CASES(typename, T, format)          \
        case typename:                          \
            PRINT_MACRO(typename, T, format);   \
            break;
        FOR_EACH_SCALAR_TYPE(CASES)
    #undef CASES
    default:
        E_RETURN(1);
        break;
    }
#undef PRINT_MACRO
    return retval;
}

static int
FmsIOAddString(FmsIOContext *ctx, const char *path, const char *value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    fprintf(ctx->fp, "%s: %s\n", path, value);
    return 0;
}

/* TODO: write the other write/read methods*/

/* Get functions. Assume that the file pointer is at the right place. */

// int
// FmsIOGetInt(FmsIOContext *ctx, const char *path, int *value)
// {
//     char *line = NULL;
//     ssize_t nread = 0;
//     size_t len = 0;
//     int retval = 0;
//     if(ctx) E_RETURN(1);
//     if(path) E_RETURN(2);
//     if(value) E_RETURN(3);
//     if((nread = getline(&line, &len, ctx->fp)) != -1)
//     {
//         size_t klen = strlen(path);
//         if(strncmp(line, path, klen) == 0)
//         {
//             if(sscanf(line+klen+1, "%d", value) != 1)
//                 retval = 4;
//         }
//         else
//             retval = 5;
//         FREE(line);
//     }
//     return retval;
// }

/**
@brief This function initializes the IO context .
@param ctx The I/O context.
*/
static void
FmsIOContextInitialize(FmsIOContext *ctx)
{
    ctx->fp = NULL;
}

/**
@brief This function initializes the functions object with function pointers that
       perform ASCII I/O.
@param ctx The functions object.
*/
static void
FmsIOFunctionsInitialize(FmsIOFunctions *obj)
{
    if(obj)
    {
        obj->open = FmsIOOpen;
        obj->close = FmsIOClose;
        obj->add_int = FmsIOAddInt;
        obj->add_int_array = FmsIOAddIntArray;
        obj->add_typed_int_array = FmsIOAddTypedIntArray;
        obj->add_float = FmsIOAddFloat;
        obj->add_double = FmsIOAddDouble;
        obj->add_scalar_array = FmsIOAddScalarArray;
        obj->add_string = FmsIOAddString;
        /*...*/
        /*obj->get_int = FmsIOGetInt; */
#ifdef MEANDERING_THOUGHTS
        obj->begin_list = FmsIOBeginList;
        obj->end_list = FmsIOEndList;
        obj->begin_list_item = FmsIOBeginListItem;
        obj->end_list_item = FmsIOEndListItem;
#endif
        /*TODO: fill in rest. */
    }
}

/*****************************************************************************
*** FMS I/O Functions for Conduit
*****************************************************************************/
#ifdef FMS_HAVE_CONDUIT
static int
FmsIOOpenConduit(FmsIOContext *ctx, const char *filename, const char *mode)
{
    if(!ctx) E_RETURN(1);
    if(!filename) E_RETURN(2);
    if(!mode) E_RETURN(3);
    ctx->filename = strdup(filename);
    ctx->root = conduit_node_create();
    ctx->writing = (strstr(mode, "w") != NULL) ? 1 : 0;

    if(!ctx->writing)
    {
        /* If we're not writing then read it all upfront into ctx->root. */
        conduit_relay_io_load(filename, ctx->protocol, NULL, ctx->root);
    }

    return 0;
}

static int
FmsIOCloseConduit(FmsIOContext *ctx)
{
    if(!ctx) E_RETURN(1);

    if(ctx->writing)
    {
        /* For conduit, we have not written anything yet. Do it all now. */
        conduit_relay_io_save(ctx->root, ctx->filename, ctx->protocol, NULL);
    }

    free(ctx->filename);
    ctx->filename = NULL;
    free(ctx->protocol);
    ctx->protocol = NULL;
    conduit_node_destroy(ctx->root);
    ctx->root = NULL;
    return 0;
}

static int
FmsIOHasPathConduit(FmsIOContext *ctx, const char *path)
{
    if(!ctx) E_RETURN(0);
    if(!path) E_RETURN(0);
    return conduit_node_has_path(ctx->root, path);
}

static int
FmsIOAddIntConduit(FmsIOContext *ctx, const char *path, FmsInt value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    /*conduit_node_set_path_int(ctx->root, path, value);*/
    conduit_node_set_path_uint64(ctx->root, path, value);
    return 0;
}

static int
FmsIOAddIntArrayConduit(FmsIOContext *ctx, const char *path, const FmsInt *values, size_t n)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);

    char *ksize = join_keys(path, "Size");
    conduit_node_set_path_uint64(ctx->root, ksize, n);
    FREE(ksize);

    char *ktype = join_keys(path, "Type");
    conduit_node_set_path_char8_str(ctx->root, ktype, FmsIntTypeNames[FMS_UINT64]);
    FREE(ktype);
    
    if(!values || n == 0)
        return 0;

    char *kvalues = join_keys(path, "Values");
    conduit_node_set_path_uint64_ptr(ctx->root, kvalues, (FmsInt *)values, n);
    FREE(kvalues);
    return 0;
}

static int
FmsIOAddTypedIntArrayConduit(FmsIOContext *ctx, const char *path, FmsIntType type, const void *values, FmsInt n)
{
    int retval = 0;
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    
    char *ksize = join_keys(path, "Size");
    conduit_node_set_path_uint64(ctx->root, ksize, n);
    FREE(ksize);

    char *ktype = join_keys(path, "Type");
    conduit_node_set_path_char8_str(ctx->root, ktype, FmsIntTypeNames[FMS_UINT64]);
    FREE(ktype);
    
    if(!values || n == 0)
        return 0;

    char *kvalues = join_keys(path, "Values");
    switch(type)
    {
    case FMS_INT8:
        conduit_node_set_path_int8_ptr(ctx->root, kvalues, (conduit_int8 *)values, n);
        break;
    case FMS_INT16:
        conduit_node_set_path_int16_ptr(ctx->root, kvalues, (conduit_int16 *)values, n);
        break;
    case FMS_INT32:
        conduit_node_set_path_int32_ptr(ctx->root, kvalues, (conduit_int32 *)values, n);
        break;
    case FMS_INT64:
        conduit_node_set_path_int64_ptr(ctx->root, kvalues, (conduit_int64 *)values, n);
        break;
    case FMS_UINT8:
        conduit_node_set_path_uint8_ptr(ctx->root, kvalues, (conduit_uint8 *)values, n);
        break;
    case FMS_UINT16:
        conduit_node_set_path_uint16_ptr(ctx->root, kvalues, (conduit_uint16 *)values, n);
        break;
    case FMS_UINT32:
        conduit_node_set_path_uint32_ptr(ctx->root, kvalues, (conduit_uint32 *)values, n);
        break;
    case FMS_UINT64:
        conduit_node_set_path_uint64_ptr(ctx->root, kvalues, (conduit_uint64 *)values, n);
        break;
    default:
        FREE(kvalues);
        E_RETURN(1);
        break;
    }
    FREE(kvalues);
    return retval;
}

static int
FmsIOAddFloatConduit(FmsIOContext *ctx, const char *path, float value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    conduit_node_set_path_float(ctx->root, path, value);
    return 0;
}

static int
FmsIOAddDoubleConduit(FmsIOContext *ctx, const char *path, double value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    conduit_node_set_path_double(ctx->root, path, value);
    return 0;
}

static int
FmsIOAddScalarArrayConduit(FmsIOContext *ctx, const char *path, FmsScalarType type, const void *data, FmsInt size)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);

    char *ksize = join_keys(path, "Size");
    conduit_node_set_path_uint64(ctx->root, ksize, size);
    FREE(ksize);

    char *ktype = join_keys(path, "Type");
    conduit_node_set_path_char8_str(ctx->root, ktype, FmsScalarTypeNames[type]);
    FREE(ktype);

    if(!data || size == 0)
        return 0;

    char *kvalues = join_keys(path, "Values");
    switch (type) {
        case FMS_FLOAT:
            conduit_node_set_path_float32_ptr(ctx->root, kvalues, (conduit_float32 *)data, size);
            break;
        case FMS_DOUBLE:
            conduit_node_set_path_float64_ptr(ctx->root, kvalues, (conduit_float64 *)data, size);
            break;
        case FMS_COMPLEX_FLOAT:
            conduit_node_set_path_float32_ptr(ctx->root, kvalues, (conduit_float32 *)data, size*2);
            break;
        case FMS_COMPLEX_DOUBLE:
            conduit_node_set_path_float64_ptr(ctx->root, kvalues, (conduit_float64 *)data, size*2);
            break;
        default:
            FREE(kvalues);
            E_RETURN(1);
            break;
    }
    FREE(kvalues);
    return 0;
}

static int
FmsIOAddStringConduit(FmsIOContext *ctx, const char *path, const char *value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    conduit_node_set_path_char8_str(ctx->root, path, value);
    /* NOT IMPLEMENTED IN CONDUIT conduit_node_set_path_external_char8_str(ctx->root, path, (char *)value);*/
    return 0;
}

static int
FmsIOGetIntConduit(FmsIOContext *ctx, const char *path, FmsInt *value)
{
    conduit_node *node = NULL;
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    if(!value) E_RETURN(3);
    if(!conduit_node_has_path(ctx->root, path))
        E_RETURN(4);

    if((node = conduit_node_fetch(ctx->root, path)) != NULL)
    {
        const conduit_datatype *dt = conduit_node_dtype(node);
        if(dt)
        {
            if(conduit_datatype_is_number(dt))
            {
                *value = (FmsInt)conduit_node_fetch_path_as_int64(ctx->root, path);
            }
            else
            {
                E_RETURN(5);
            }
        }
    }
    return 0;
}

/*int (*get_typed_int_array)(FmsIOContext *ctx, const char *path, FmsIntType *type, void **values, FmsInt *n);*/
static int
FmsIOGetTypedIntArrayConduit(FmsIOContext *ctx, const char *path, FmsIntType *type, void **values, FmsInt *n)
{
    conduit_node *node = NULL;
    char *type_str = NULL;
    void *conduit_data_ptr = NULL;
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    if(!type) E_RETURN(3);
    if(!values) E_RETURN(4);
    if(!n) E_RETURN(5);

    if(!conduit_node_has_path(ctx->root, path))
        E_RETURN(6);

    if((node = conduit_node_fetch(ctx->root, path)) == NULL)
        E_RETURN(7);

    /* Arrays are written with subnodes Size, Type, and Values */
    if(!conduit_node_has_child(node, "Size"))
        E_RETURN(8);
    *n = (FmsInt)conduit_node_fetch_path_as_int64(node, "Size");

    if(!conduit_node_has_child(node, "Type"))
        E_RETURN(9);
    
    if((type_str = conduit_node_fetch_path_as_char8_str(node, "Type")) == NULL)
        E_RETURN(10);

    if(FmsGetIntTypeFromName(type_str, type))
        E_RETURN(11);
    
    if(*n == 0)
    {
        *values = NULL;
        return 0;
    }

    if(!conduit_node_has_child(node, "Values"))
        E_RETURN(12);

    conduit_data_ptr = conduit_node_element_ptr(conduit_node_fetch(node, "Values"), 0);

    /* TODO fix this - can't do a switch over our types because conduit is probably using a signed pointer.
        Maybe use conduit_datatype_is_signed()
        or do if/else if over each conduit_datatype for the correct one  */
    switch(*type)
    {
    #define COPY_DATA(typename, T, format) \
        case typename: \
        { \
            void *dest = malloc(sizeof(T) * *n); \
            memcpy(dest, conduit_data_ptr, sizeof(T) * *n); \
            *values = dest; \
            break; \
        }
        FOR_EACH_INT_TYPE(COPY_DATA)
    #undef COPY_DATA
        default:
            E_RETURN(13);
    }

    return 0;
}

/* int (*get_scalar_array)(FmsIOContext *ctx, const char *path, FmsScalarType *type, void **values, FmsInt *n); */
static int
FmsIOGetScalarArrayConduit(FmsIOContext *ctx, const char *path, FmsScalarType *type, void **values, FmsInt *n)
{
    conduit_node *node = NULL;
    char *type_str = NULL;
    void *conduit_data_ptr = NULL;
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    if(!type) E_RETURN(3);
    if(!values) E_RETURN(4);
    if(!n) E_RETURN(5);

    /* Arrays are written with subnodes Size, Type, and Values */
    if((node = conduit_node_fetch(ctx->root, path)) == NULL)
        E_RETURN(6);

    if((type_str = conduit_node_fetch_path_as_char8_str(node, "Type")) == NULL)
        E_RETURN(7);

    if(FmsGetScalarTypeFromName(type_str, type))
        E_RETURN(8);

    
    if(!conduit_node_has_child(node, "Size"))
        E_RETURN(9);

    *n = (FmsInt)conduit_node_fetch_path_as_int64(node, "Size");

    if(n == 0)
    {
        *values = NULL;
        return 0;
    }

    if(!conduit_node_has_child(node, "Values"))
        E_RETURN(10);

    conduit_data_ptr = conduit_node_element_ptr(conduit_node_fetch(node, "Values"), 0);

    /* Need to copy the data out of conduit - right? */
    switch(*type)
    {
    #define COPY_DATA(typename, T, format) \
        case typename: \
        { \
            void *dest = malloc(sizeof(T) * *n); \
            memcpy(dest, conduit_data_ptr, sizeof(T) * *n); \
            *values = dest; \
            break; \
        }
        FOR_EACH_SCALAR_TYPE(COPY_DATA)
    #undef COPY_DATA
        default:
            E_RETURN(11);
    }

    return 0;
}

/* int (*get_string)(FmsIOContext *ctx, const char *path, char **value); */
static int
FmsIOGetStringConduit(FmsIOContext *ctx, const char *path, const char **value)
{
    conduit_node *node = NULL;
    const conduit_datatype *dt = NULL;
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    if(!value) E_RETURN(3);

    if(!conduit_node_has_path(ctx->root, path))
        E_RETURN(4);

    if((node = conduit_node_fetch(ctx->root, path)) == NULL)
        E_RETURN(5);

    dt = conduit_node_dtype(node);

    if(!conduit_datatype_is_char8_str(dt))
        E_RETURN(6);

    *value = conduit_node_fetch_path_as_char8_str(ctx->root, path);
    
    return 0;
}

/*  int (*get_string_array)(FmsIOContext *ctx, const char *path, char ***value, FmsInt *n); */
static int
FmsIOGetStringArrayConduit(FmsIOContext *ctx, const char *path, const char ***values, FmsInt *n)
{
    /* TODO: */
    E_RETURN(1);
}

#ifdef MEANDERING_THOUGHTS
/**
@brief Push the conduit node to the top of the stack, making it current.
@param ctx  The context.
@param node The conduit node to push onto the stack.
*/
static void
FMSIOContextPush(FmsIOContext *ctx, conduit_node *node)
{
    if(ctx->stack_top+1 < 100)
    {
        ctx->stack_top++;
        ctx->stack[ctx->stack_top] = node;
    }
}

/**
@brief Pop the top Conduit node off the stack, making the one "under" it current.
@param ctx  The context.
*/
static void
FMSIOContextPop(FmsIOContext *ctx)
{
    if(ctx->stack_top > -1)
    {
        ctx->stack[ctx->stack_top] = NULL;
        ctx->stack_top--;
    }
}

/**
@brief Get the current Conduit node on the top of the stack. Keys are made relative 
       to this node. This means that except for the root node, we should use relative
       paths.
@param ctx  The context.
@return Return the current Conduit node in the context.
*/
static conduit_node *
FmsIOContextCurrent(FmsIOContext *ctx)
{
    conduit_node *current = ctx->root;
    if(ctx->stack_top > -1)
        current = ctx->stack[ctx->stack_top];
    return current;
}

static int
FmsIOBeginListConduit(FmsIOContext *ctx, const char *path, size_t n)
{
    conduit_node *node = NULL;
    if(ctx) E_RETURN(1);
    if(path) E_RETURN(2);
    if(value) E_RETURN(3);

    /* Get a new node relative to the current, creating it if necessary. */
    if((node = conduit_node_fetch(FMSIOContextCurrent(ctx), path)) != NULL)
    {
        /* Make the new node the current node. */
        FmsIOContextPush(ctx, node);
        conduit_node *conduit_node_append(conduit_node *cnode);
    }
}

static int
FmsIOEndListConduit(FmsIOContext *ctx)
{
    if(ctx) E_RETURN(1);
    FmsIOContextPop(ctx);
    return 0;
}

static int
FmsIOBeginListItemConduit(FmsIOContext *ctx)
{
    conduit_node *node = NULL;
    if(ctx) E_RETURN(1);
    if(path) E_RETURN(2);
    if(value) E_RETURN(3);

    /* Get a new node relative to the current*/
    node = conduit_node_append(conduit_node_create());
    FMSIOContextPush(ctx, node);
}

static int
FmsIOEndListItemConduit(FmsIOContext *ctx)
{
    if(ctx) E_RETURN(1);
    FmsIOContextPop(ctx);
    return 0;
}
#endif

/* TODO: write the other write/read methods*/




/**
@brief This function initializes the IO context with function pointers that
       perform ASCII I/O.
@param ctx The I/O context.
*/
static void
FmsIOContextInitializeConduit(FmsIOContext *ctx, const char *protocol)
{
    if(ctx)
    {
        int i;
#ifdef MEANDERING_THOUGHTS
        for(i = 0; i < 100; ++i)
            ctx->stack[i] = NULL;
        ctx->stack_top = -1;
#endif

        ctx->filename = NULL;
        ctx->root = NULL;
        ctx->writing = 0;
        /* Store the protocol this way since we do not pass it to the open function */
        ctx->protocol = protocol ? strdup(protocol) : strdup("json");
    }
}

/**
@brief This function initializes the functions object with function pointers that
       perform ASCII I/O.
@param ctx The functions object.
*/
static void
FmsIOFunctionsInitializeConduit(FmsIOFunctions *obj)
{
    if(obj)
    {
        obj->open = FmsIOOpenConduit;
        obj->close = FmsIOCloseConduit;
        obj->add_int = FmsIOAddIntConduit;
        obj->add_int_array = FmsIOAddIntArrayConduit;
        obj->add_typed_int_array = FmsIOAddTypedIntArrayConduit;
        obj->add_float = FmsIOAddFloatConduit;
        obj->add_double = FmsIOAddDoubleConduit;
        obj->add_string = FmsIOAddStringConduit;
        obj->add_scalar_array = FmsIOAddScalarArrayConduit;
        /*...*/
        obj->has_path = FmsIOHasPathConduit;
        obj->get_int = FmsIOGetIntConduit;
        obj->get_typed_int_array = FmsIOGetTypedIntArrayConduit;
        obj->get_scalar_array = FmsIOGetScalarArrayConduit;
        obj->get_string = FmsIOGetStringConduit;
        obj->get_string_array = FmsIOGetStringArrayConduit;
#ifdef MEANDERING_THOUGHTS
        obj->begin_list = FmsIOBeginListConduit;
        obj->end_list = FmsIOEndListConduit;
        obj->begin_list_item = FmsIOBeginListItemConduit;
        obj->end_list_item = FmsIOEndListItemConduit;
#endif
        /*TODO: fill in rest. */
    }
}


#endif

/*****************************************************************************
*** FMS I/O Building Block Functions
*****************************************************************************/

/**
@brief Write FmsMetaData to the output I/O context.
@param ctx The context
@param io  The I/O functions that operate on the context.
@param key The key that identifies the data in the output.
@param mdata  The FMS metadata.
@return 0 on success, non-zero otherwise.
*/
static int
FmsIOWriteFmsMetaData(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsMetaData mdata)
{
    /** TODO: write me */
    int err = 0;
    FmsMetaDataType type;

    /* It says it returns 1 if mdata is empty, but then just checks if its null (if im not mistaken).
        We already null checked before we got here so I'm not sure if we need to check for 1. */        
    int r = FmsMetaDataGetType(mdata, &type);
    if(r > 1)
        E_RETURN(1);
    else if(r == 1)
        return 0;

    char *kmdtype = join_keys(key, "MetaDataType");
    (*io->add_string)(ctx, kmdtype, FmsMetaDataTypeNames[type]);
    FREE(kmdtype);

    switch(type)
    {
        case FMS_INTEGER:
        {
            char *kname = NULL, *kdata = NULL;
            const char *mdata_name = NULL;
            FmsIntType int_type;
            FmsInt size;
            const void *data = NULL;

            if(FmsMetaDataGetIntegers(mdata, &mdata_name, &int_type, &size, &data))
                E_RETURN(3);
            
            kname = join_keys(key, "Name");
            err = (*io->add_string)(ctx, kname, mdata_name);
            FREE(kname);

            kdata = join_keys(key, "Data");
            err = (*io->add_typed_int_array)(ctx, kdata, int_type, data, size);
            FREE(kdata);                
            if(err)
                E_RETURN(4);
            break;
        }
        case FMS_SCALAR:
        {
            char *kname = NULL, *kdata = NULL;
            const char *mdata_name = NULL;
            FmsScalarType scalar_type;
            FmsInt size;
            const void *data = NULL;

            if(FmsMetaDataGetScalars(mdata, &mdata_name, &scalar_type, &size, &data))
                E_RETURN(5);
            
            kname = join_keys(key, "Name");
            err = (*io->add_string)(ctx, kname, mdata_name);
            FREE(kname);

            kdata = join_keys(key, "Data");
            err = (*io->add_scalar_array)(ctx, kdata, scalar_type, data, size);
            FREE(kdata);                
            if(err)
                E_RETURN(6);
            break;
        }
        case FMS_STRING:
        {
            char *kname = NULL, *kdata = NULL;
            const char *mdata_name = NULL;
            const char *data = NULL;

            if(FmsMetaDataGetString(mdata, &mdata_name, &data))
                E_RETURN(7);
            
            kname = join_keys(key, "Name");
            err = (*io->add_string)(ctx, kname, mdata_name);
            FREE(kname);

            kdata = join_keys(key, "Data");
            err = (*io->add_string)(ctx, kdata, data);
            FREE(kdata);       
            if(err)
                E_RETURN(8);
            break;
        }
        case FMS_META_DATA:
        {
            char *kname = NULL, *kdata = NULL, *ksize = NULL;
            const char *mdata_name = NULL;
            FmsInt i, size = 0;
            FmsMetaData *data = NULL;

            if(FmsMetaDataGetMetaData(mdata, &mdata_name, &size, &data))
                E_RETURN(9);

            if(!data)
                E_RETURN(10);

            kname = join_keys(key, "Name");
            err = (*io->add_string)(ctx, kname, mdata_name);
            FREE(kname);
            
            ksize = join_keys(key, "Size");
            err = (*io->add_int)(ctx, ksize, size);
            FREE(ksize);

            kdata = join_keys(key, "Data");
            // Recursive?
            for(i = 0; i < size; i++)
            {
                char temp[20], *tk = NULL;
                sprintf(temp, "%d", (int)i);
                tk = join_keys(kdata, temp);
                err = FmsIOWriteFmsMetaData(ctx, io, tk, data[i]);
            }
            FREE(kdata);

            if(err)
                E_RETURN(11);
            break;
        }
        default:
            E_RETURN(12);
            break;
    }        
    return 0;         
}

static int
FmsIOReadFmsMetaData(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, 
    FmsIOMetaDataInfo *mdinfo)
{
    int err = 0;
    if(!mdinfo) E_RETURN(1);

    {
        char *kmdtype = join_keys(key, "MetaDataType");
        const char *str = NULL;
        err = (*io->get_string)(ctx, kmdtype, &str);
        FREE(kmdtype);
        if(err)
            E_RETURN(2);
        if(FmsGetMetaDataTypeFromName(str, &mdinfo->mdtype))
            E_RETURN(3);
    }

    // Get name
    {
        char *kname = join_keys(key, "Name");
        err = (*io->get_string)(ctx, kname, &mdinfo->name);
        FREE(kname);
        if(err)
            E_RETURN(4);
    }

    switch (mdinfo->mdtype)
    {
        case FMS_INTEGER:
        {
            // Get DataArray
            char *kdata = join_keys(key, "Data");
            err = (*io->get_typed_int_array)(ctx, kdata, &mdinfo->i_type, &mdinfo->data, &mdinfo->size);
            FREE(kdata);
            if(err)
                E_RETURN(5);
            break;
        }
        case FMS_SCALAR:
        {
            // Get data array
            char *kdata = join_keys(key, "Data");
            err = (*io->get_scalar_array)(ctx, kdata, &mdinfo->s_type, &mdinfo->data, &mdinfo->size);
            FREE(kdata);
            if(err)
                E_RETURN(6);
            break;
        }
        case FMS_STRING:
        {
            // Get data array
            char *kdata = join_keys(key, "Data");
            err = (*io->get_string_array)(ctx, kdata, (const char***)&mdinfo->data, &mdinfo->size);
            FREE(kdata);
            if(err)
                E_RETURN(7);
            break;
        }
        case FMS_META_DATA:
        {
            // Get size
            FmsInt i;
            FmsIOMetaDataInfo *mds = NULL;
            char *ksize = join_keys(key, "Size"), *kdata = NULL;
            err = (*io->get_int)(ctx, ksize, &mdinfo->size);
            FREE(ksize);
            if(err)
                E_RETURN(8);
            mdinfo->data = calloc(sizeof(FmsIOMetaDataInfo), mdinfo->size);
            mds = (FmsIOMetaDataInfo *)mdinfo->data;
            kdata = join_keys(key, "Data");
            for(i = 0; i < mdinfo->size; i++)
            {
                char temp[21], *ki = NULL;
                sprintf(temp, FMS_LU, i);
                ki = join_keys(kdata, temp);
                err |= FmsIOReadFmsMetaData(ctx, io, ki, &mds[i]);
            }
            FREE(kdata);
            if(err)
                E_RETURN(9);
            break;
        }
        default:
            break;
    }
    return 0;
}

static int
FmsIOWriteFmsDomain(FmsIOContext *ctx, FmsIOFunctions *io, const char *key,
    FmsDomain dom)
{
    FmsInt dim = 0;
    if(FmsDomainGetDimension(dom, &dim))
        E_RETURN(1);
    char *kdim = join_keys(key, "Dimension");
    (io->add_int)(ctx, kdim, dim);
    FREE(kdim);

    FmsInt nverts = 0;
    if(FmsDomainGetNumVertices(dom, &nverts))
        E_RETURN(2);
    char *knverts = join_keys(key, "NumVertices");
    (io->add_int)(ctx, knverts, nverts);
    FREE(knverts);

    // Skip verticies because we already know what we need to know
    char *kentities = join_keys(key, "Entities");
    FmsInt i, entries = 0;
    for(i = 1; i < FMS_NUM_ENTITY_TYPES; i++)
    {
        FmsInt nents = 0;
        FmsDomainGetNumEntities(dom, (FmsEntityType)i, &nents);
        if(nents)
        {
            char temp[8];
            sprintf(temp, "%d", (int)entries++);
            char *kentry = join_keys(kentities, temp);
            const char *enttype = FmsEntityTypeNames[i];
            char *kenttype = join_keys(kentry, "EntityType");
            (*io->add_string)(ctx, kenttype, enttype);
            FREE(kenttype);

            char *kne = join_keys(kentry, "NumEntities");
            (*io->add_int)(ctx, kne, nents);
            FREE(kne);

            FmsIntType idtype;
            const void *ents = NULL;
            FmsInt ne = 0;
            if(FmsDomainGetAllEntities(dom, (FmsEntityType)i, 
                &idtype, &ents, &ne))
            {
                FREE(kentry); E_RETURN(3);
            }

            if(ents == NULL)
            {
                FREE(kentry); E_RETURN(4);
            }

            // TODO: Remove this check
            if(nents != ne)
            {
                FREE(kentry); E_RETURN(5);
            }

            const FmsInt N = FmsEntityNumSides[i] * nents;
            (*io->add_typed_int_array)(ctx, kentry, idtype, ents, N);
            FREE(kentry);
        }
    }
    FREE(kentities);
    // Q: Do we even have to store the orientations? Fms converts everything to standard Fms orientation in memeory.
    return 0;
}

static int
FmsIOReadFmsDomain(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsIODomainInfo *dinfo) {
    int err = 0;
    if(!dinfo)
        E_RETURN(1);

    // Get dimension
    {
        char *kdom = join_keys(key, "Dimension");
        err = (*io->get_int)(ctx, kdom, &dinfo->dim);
        FREE(kdom);
        if(err)
            E_RETURN(2);
    }

    // Get Num verticies
    {
        char *knverts = join_keys(key, "NumVertices");
        err = (*io->get_int)(ctx, knverts, &dinfo->nverts);
        FREE(knverts);
        if(err)
            E_RETURN(3);
    }

    if(dinfo->dim)
    {
        char *kentities = join_keys(key, "Entities");
        FmsInt i = 0;
        dinfo->entities = calloc(sizeof(FmsIOEntityInfo), dinfo->dim);
        for(i = 0; i < dinfo->dim; i++)
        {
            char temp[21], *ke = NULL;
            sprintf(temp, FMS_LU, i);
            ke = join_keys(kentities, temp);

            // Get EntityType
            {
                char *kent_type = join_keys(ke, "EntityType");
                const char *ent_name = NULL;
                err = (*io->get_string)(ctx, kent_type, &ent_name);
                FREE(kent_type);
                if(err)
                    E_RETURN(4);
                if(FmsGetEntityTypeFromName(ent_name, &dinfo->entities[i].ent_type))
                    E_RETURN(5);
            }

            // Get NumEnts
            {
                char *knents = join_keys(ke, "NumEntities");
                err = (*io->get_int)(ctx, knents, &dinfo->entities[i].nents);
                FREE(knents);
                if(err)
                    E_RETURN(6);
            }

            // Get Entity array
            err = (*io->get_typed_int_array)(ctx, ke, &dinfo->entities[i].type, 
                &dinfo->entities[i].values, &dinfo->entities[i].size);
            FREE(ke);
            if(err)
                E_RETURN(7);
        }
    }
    return 0;
}

static int
FmsIOWriteDomains(FmsIOContext *ctx, FmsIOFunctions *io, const char *key,
    FmsMesh mesh, FmsInt di)
{
    const char *domname = NULL;
    FmsInt ndoms = 0;
    FmsDomain *doms = NULL;
    if(FmsMeshGetDomains(mesh, di, &domname, &ndoms, &doms))
        E_RETURN(1);

    if(!domname)
        E_RETURN(2);
    char *kdname = join_keys(key, "Name");
    (*io->add_string)(ctx, kdname, domname);
    FREE(kdname);

    char *knd = join_keys(key, "NumDomains");
    (*io->add_int)(ctx, knd, ndoms);
    FREE(knd);

    char *kdoms = join_keys(key, "Domains");
    FmsInt i;
    for(i = 0; i < ndoms; i++)
    {
        char tmp[20], *kd = NULL;
        sprintf(tmp, "%d", (int)i);
        kd = join_keys(kdoms, tmp);
        FmsIOWriteFmsDomain(ctx, io, kd, doms[i]);
        FREE(kd);
    }
    FREE(kdoms);
    return 0;
}

static int
FmsIOReadFmsDomainName(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsIODomainNameInfo *dinfo) {
    int err = 0;
    if(!dinfo)
        E_RETURN(1);

    // Get domain name
    {
        char *kdname = join_keys(key, "Name");
        err = (*io->get_string)(ctx, kdname, &dinfo->name);
        FREE(kdname);
        if(err)
            E_RETURN(2);
    }

    // Get number of domains
    {
        char *kndoms = join_keys(key, "NumDomains");
        err = (*io->get_int)(ctx, kndoms, &dinfo->ndomains);
        FREE(kndoms);
        if(err)
            E_RETURN(3);
    }

    // Read the domains (if we have any)
    if(dinfo->ndomains)
    {
        char *kdoms = join_keys(key, "Domains");
        FmsInt i = 0;
        dinfo->domains = calloc(sizeof(FmsIODomainInfo), dinfo->ndomains);
        for(i = 0; i < dinfo->ndomains; i++)
        {
            char temp[21], *kd = NULL;
            sprintf(temp, FMS_LU, i);
            kd = join_keys(kdoms, temp);
            err |= FmsIOReadFmsDomain(ctx, io, kd, &dinfo->domains[i]);
            FREE(kd);
        }
        if(err)
            E_RETURN(4);
    }
    return 0;
}

/**
@brief Write FmsMesh to the output I/O context.
@param ctx The context
@param io  The I/O functions that operate on the context.
@param key The key that identifies the data in the output.
@param comp  The FMS component.
@return 0 on success, non-zero otherwise.
*/
static int
FmsIOWriteComponent(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsComponent comp)
{
    const char *comp_name = NULL;
    if(FmsComponentGetName(comp, &comp_name))
        E_RETURN(1);
    char *kcomp_name = join_keys(key, "Name");
    (*io->add_string)(ctx, kcomp_name, comp_name);
    FREE(kcomp_name);

    FmsInt dim = 0;
    if(FmsComponentGetDimension(comp, &dim))
        E_RETURN(2);
    char *kdim = join_keys(key, "Dimension");
    (*io->add_int)(ctx, kdim, dim);
    FREE(kdim);

    FmsInt nents = 0;
    if(FmsComponentGetNumEntities(comp, &nents))
        E_RETURN(3);
    char *knents = join_keys(key, "NumEntities");
    (*io->add_int)(ctx, knents, nents);
    FREE(knents);

    FmsField coords = NULL;
    if(FmsComponentGetCoordinates(comp, &coords))
        E_RETURN(4);
    // Coords can be null
    if(coords)
    {
        const char *coords_name = NULL;
        FmsFieldGetName(coords, &coords_name);
        char *kcoords_name = join_keys(key, "Coordinates");
        (*io->add_string)(ctx, kcoords_name, coords_name);
        FREE(kcoords_name);
    }

    FmsInt nparts = 0;
    if(FmsComponentGetNumParts(comp, &nparts))
        E_RETURN(5);
    char *knparts = join_keys(key, "NumParts");
    (*io->add_int)(ctx, knparts, nparts);
    FREE(knparts);

    char *kparts = join_keys(key, "Parts");
    FmsInt i;
    for(i = 0; i < nparts; i++)
    {
        char temp[20], *kpart;
        sprintf(temp, "%d", (int)i);
        kpart = join_keys(kparts, temp);
        FmsDomain dom;
        if(FmsComponentGetPart(comp, i, 0, &dom, NULL, NULL, NULL, NULL))
            E_RETURN(6);

        if(dom)
        {
            const char *domain_name = NULL;
            FmsInt did = 0;
            if(FmsDomainGetName(dom, &domain_name, &did))
                E_RETURN(7);
            
            char *kdomain_name = join_keys(kpart, "DomainName");
            (*io->add_string)(ctx, kdomain_name, domain_name);
            FREE(kdomain_name);

            char *kdid = join_keys(kpart, "DomainID");
            (*io->add_int)(ctx, kdid, did);
            FREE(kdid);
        }
        
        char found = 0;
        FmsInt et;
        for(et = 0; et < FMS_NUM_ENTITY_TYPES; et++)
        {
            FmsIntType ent_id_type;
            const void *ents;
            FmsInt nents = 0;
            // Q: Doing the for loop like this should pickup all the local indexing too?
            //  FmsComponentGetPart and FmsComponentGetPartSubEntities index into the same arrays.
            if(FmsComponentGetPart(comp, i, (FmsEntityType)et, 
                NULL, &ent_id_type, &ents, NULL, &nents))
                E_RETURN(8);
            // NULL means use all entities in that domain
            if(ents && nents)
            {
                found = 1;
                const char *enttype = FmsEntityTypeNames[et];
                if(FmsEntityDim[et] == dim)
                {
                    char *kenttype = join_keys(kpart, "PartEntityType");
                    (*io->add_string)(ctx, kenttype, enttype);
                    FREE(kenttype);
                }
                char *k = join_keys(kpart, enttype);
                
                char *knents = join_keys(k, "NumEntities");
                (*io->add_int)(ctx, knents, nents);
                FREE(knents);

                const FmsInt N = FmsEntityNumSides[et] * nents;
                (*io->add_typed_int_array)(ctx, k, ent_id_type, ents, N);
                FREE(k);
            }
        }
        // "Typically, the whole is represented by a special component of all elements (3-entities) on all domains."
        if(!found)
        {
            char *kfull_domain = join_keys(kpart, "FullDomain");
            (*io->add_string)(ctx, kfull_domain, "Yes");
            FREE(kfull_domain);            
        }
        FREE(kpart);
    }
    FREE(kparts);

    const FmsInt *rel_comps = NULL;
    FmsInt nrelcomps = 0;
    if(FmsComponentGetRelations(comp, &rel_comps, &nrelcomps))
        E_RETURN(9);
    

    char *krelations = join_keys(key, "Relations");
    (*io->add_int_array)(ctx, krelations, rel_comps, nrelcomps);
    FREE(krelations);

    return 0;
}

static int
FmsIOReadFmsComponent(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsIOComponentInfo *comp_info) {
    int err = 0;
    if(!comp_info) E_RETURN(1);

    // Get Component name
    {
        char *kcomp_name = join_keys(key, "Name");
        err = (*io->get_string)(ctx, kcomp_name, &comp_info->name);
        FREE(kcomp_name);
        if(err)
            E_RETURN(2);
    }

    // Get dim
    {
        char *kcomp_dim = join_keys(key, "Dimension");
        err = (*io->get_int)(ctx, kcomp_dim, &comp_info->dim);
        FREE(kcomp_dim);
        if(err)
            E_RETURN(3);
    }

    {
        char *knents = join_keys(key, "NumEntities");
        err = (*io->get_int)(ctx, knents, &comp_info->nents);
        FREE(knents);
        if(err)
            E_RETURN(4);
    }

    /// Not all components have Coords
    {
        char *kcoords = join_keys(key, "Coordinates");
        if((*io->has_path)(ctx, kcoords))
            err = (*io->get_string)(ctx, kcoords, &comp_info->coordinates_name);
        FREE(kcoords);
        if(err)
            E_RETURN(5);
    }

    {
        char *knparts = join_keys(key, "NumParts");
        err = (*io->get_int)(ctx, knparts, &comp_info->nparts);
        FREE(knparts);
        if(err)
            E_RETURN(6);
    }

    // If we have parts get their info
    if(comp_info->nparts)
    {
        char *kparts = join_keys(key, "Parts");
        FmsInt i, et;
        comp_info->parts = calloc(sizeof(FmsIOComponentPartInfo), comp_info->nparts);
        for(i = 0; i < comp_info->nparts; i++)
        {
            char temp[21], *kpart = NULL;
            sprintf(temp, FMS_LU, i);
            kpart = join_keys(kparts, temp);
            
            // Get the DomainName for this part
            {
                char *kdom_name = join_keys(kpart, "DomainName");
                err = (*io->get_string)(ctx, kdom_name, &comp_info->parts[i].domain_name);
                FREE(kdom_name);
                if(err)
                    E_RETURN(7);
            }

            // Get the domain id for this part
            {
                char *kdid = join_keys(kpart, "DomainID");
                err = (*io->get_int)(ctx, kdid, &comp_info->parts[i].domain_id);
                FREE(kdid);
                if(err)
                    E_RETURN(8);
            }

            // Check to see if this is a FullDomain part, this part will be complete if true
            {
                char *kfulldomain = join_keys(kpart, "FullDomain");
                if((*io->has_path)(ctx, kfulldomain))
                {
                    const char *isfulldomain = NULL;
                    err = (*io->get_string)(ctx, kfulldomain, &isfulldomain);
                    if(err)
                        E_RETURN(9);
                    // Special case, this part is now complete
                    if(strcmp("Yes", isfulldomain) == 0)
                    {
                        comp_info->parts[i].full_domain = 1;
                        FREE(kfulldomain); FREE(kpart);
                        continue;
                    }
                }
                comp_info->parts[i].full_domain = 0;
                FREE(kfulldomain);
            }

            // If this was not a full domain then we have to figure out the entity & subentities
            {
                char *kentity = join_keys(kpart, "PartEntityType");
                const char *str = NULL;
                err = (io->get_string)(ctx, kentity, &str);
                FREE(kentity);
                if(err)
                    E_RETURN(10);
                if(FmsGetEntityTypeFromName(str, &comp_info->parts[i].part_ent_type))
                    E_RETURN(11);
            }

            // Entity type is encoded into the key for these
            for(et = 0; et < FMS_NUM_ENTITY_TYPES; et++)
            {
                char *k = join_keys(kpart, FmsEntityTypeNames[et]);
                if((*io->has_path)(ctx, k))
                {
                    // Now we have to populate the entitiy info
                    comp_info->parts[i].entities[et].ent_type = et;
                
                    // Get nents
                    {
                        char *knents = join_keys(k, "NumEntities");
                        err = (*io->get_int)(ctx, knents, &comp_info->parts[i].entities[et].nents);
                        FREE(knents);
                    }

                    err |= (*io->get_typed_int_array)(ctx, k, &comp_info->parts[i].entities[et].type,
                        &comp_info->parts[i].entities[et].values, &comp_info->parts[i].entities[et].size);
                }
                FREE(k);
                if(err)
                    E_RETURN(12);
            }
            FREE(kpart);
        }
        FREE(kparts);
    }

    // Finally get the relations
    {
        char *krelcomps = join_keys(key, "Relations");
        err = (*io->get_typed_int_array)(ctx, krelcomps, &comp_info->relation_type,
            &comp_info->relation_values, &comp_info->relation_size);
        FREE(krelcomps);
        if(err)
            E_RETURN(13);
    }
    return 0;
}

/**
@brief Write FmsMesh to the output I/O context.
@param ctx The context
@param io  The I/O functions that operate on the context.
@param key The key that identifies the data in the output.
@param tag  The FMS tag.
@return 0 on success, non-zero otherwise.
*/
static int
FmsIOWriteTag(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsTag tag)
{
    const char *tag_name = NULL;
    if(FmsTagGetName(tag, &tag_name))
        E_RETURN(1);
    char *ktag_name = join_keys(key, "TagName");
    (*io->add_string)(ctx, ktag_name, tag_name);
    FREE(ktag_name);

    FmsComponent comp;
    if(FmsTagGetComponent(tag, &comp))
        E_RETURN(2);
    const char *comp_name = NULL;
    if(FmsComponentGetName(comp, &comp_name))
        E_RETURN(3);
    char *kcomp = join_keys(key, "Component");
    (*io->add_string)(ctx, kcomp, comp_name);
    FREE(kcomp);

    FmsIntType tag_type;
    const void *ent_tags = NULL;
    FmsInt ntags = 0;
    if(FmsTagGet(tag, &tag_type, &ent_tags, &ntags))
        E_RETURN(4);
    
    (*io->add_typed_int_array)(ctx, key, tag_type, ent_tags, ntags);

    const void *tagvalues = NULL;
    const char * const *descriptions = NULL;
    if(FmsTagGetDescriptions(tag, NULL, &tagvalues, &descriptions, &ntags))
        E_RETURN(5);
    
    char *kdescriptions = join_keys(key, "Descriptions");

    char *kdescriptions_size = join_keys(kdescriptions, "Size");
    (*io->add_int)(ctx, kdescriptions_size, ntags);
    FREE(kdescriptions_size);

    FmsInt i;
    for(i = 0; i < ntags; i++)
    {
        char temp[21], *kd;
        sprintf(temp, FMS_LU, i);
        kd = join_keys(kdescriptions, temp);
        
        {
            char *kv = join_keys(kd, "Value");
            // TODO: Fix add_int to support signed ints
            switch(tag_type)
            {
            #define CASES(typename, T, format) \
                case typename: \
                { \
                    T *vs = (T*)tagvalues; \
                    (*io->add_int)(ctx, kv, (FmsInt)vs[i]); \
                    break; \
                }
                FOR_EACH_INT_TYPE(CASES)
            #undef CASES
                default:
                    E_RETURN(6);
                    break;
            }
            FREE(kv);
        }
        
        {
            char *kdstring = join_keys(kd, "Description");
            (*io->add_string)(ctx, kdstring, descriptions[i]);
            FREE(kdstring);
        }
    }
    FREE(kdescriptions);
    return 0;
}

static int
FmsIOReadFmsTag(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsIOTagInfo *tinfo) {
    int err = 0;
    if(!tinfo)
        E_RETURN(1);

    {
        char *ktag_name = join_keys(key, "TagName");
        err = (*io->get_string)(ctx, ktag_name, &tinfo->name);
        FREE(ktag_name);
        if(err)
            E_RETURN(2);
    }

    {
        char *kcomp_name = join_keys(key, "Component");
        err = (*io->get_string)(ctx, kcomp_name, &tinfo->comp_name);
        FREE(kcomp_name);
        if(err)
            E_RETURN(3);
    }

    err = (*io->get_typed_int_array)(ctx, key, &tinfo->type, &tinfo->values, &tinfo->size);
    if(err)
        E_RETURN(4);

    char *kdescriptions = join_keys(key, "Descriptions");
    if((*io->has_path)(ctx, kdescriptions))
    {
        FmsInt i;
        
        {
            char *kdesc_size = join_keys(kdescriptions, "Size");
            err = (*io->get_int)(ctx, kdesc_size, &tinfo->descr_size);
            FREE(kdesc_size);
        }

        tinfo->tag_values = calloc(sizeof(FmsInt), tinfo->descr_size);
        tinfo->descriptions = calloc(sizeof(char*), tinfo->descr_size);
        for(i = 0; i < tinfo->descr_size; i++)
        {
            char temp[21], *ki = NULL, *kv = NULL, *kd = NULL;
            sprintf(temp, FMS_LU, i);
            ki = join_keys(kdescriptions, temp);

            // TODO: Support unsigned
            kv = join_keys(ki, "Value");
            err |= (*io->get_int)(ctx, kv, &tinfo->tag_values[i]);
            FREE(kv);

            kd = join_keys(ki, "Description");
            err |= (*io->get_string)(ctx, kd, &tinfo->descriptions[i]);
            FREE(kd);
            FREE(ki);
        }
    }
    FREE(kdescriptions);
    if(err)
        E_RETURN(5);
    
    return 0;
}


/**
@brief Write FmsMesh to the output I/O context.
@param ctx The context
@param io  The I/O functions that operate on the context.
@param key The key that identifies the data in the output.
@param mesh  The FMS mesh.
@return 0 on success, non-zero otherwise.
*/
static int
FmsIOWriteFmsMesh(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsMesh mesh)
{
    /** IDEA: Come back and just reuse one char (stack?) buffer for all these keys? */
    FmsInt pinfo[2];
    if(FmsMeshGetPartitionId(mesh, &pinfo[0], &pinfo[1]))
        E_RETURN(1);
    
    char *kpinfo = join_keys(key, "PartitionInfo");
    (*io->add_int_array)(ctx, kpinfo, pinfo, 2);
    FREE(kpinfo);

    FmsInt ndnames = 0;
    if(FmsMeshGetNumDomainNames(mesh, &ndnames))
        E_RETURN(2);

    char *kndnames = join_keys(key, "NumDomainNames");
    (*io->add_int)(ctx, kndnames, ndnames);
    FREE(kndnames);

    FmsInt ncomps = 0;
    if(FmsMeshGetNumComponents(mesh, &ncomps))
        E_RETURN(3);
    char *kncomps = join_keys(key, "NumComponents");
    (*io->add_int)(ctx, kncomps, ncomps);
    FREE(kncomps);

    FmsInt ntags = 0;
    if(FmsMeshGetNumTags(mesh, &ntags))
        E_RETURN(4);
    char *kntags = join_keys(key, "NumTags");
    (*io->add_int)(ctx, kntags, ntags);
    FREE(kntags);

    char *kdom_names = join_keys(key, "DomainNames");
    FmsInt i;
    for(i = 0; i < ndnames; i++)
    {
        char tmp[20], *dk = NULL;
        sprintf(tmp, "%d", (int)i);
        dk = join_keys(kdom_names, tmp);
        FmsIOWriteDomains(ctx, io, dk, mesh, i);
        FREE(dk);
    }
    FREE(kdom_names);

    char *kcomps = join_keys(key, "Components");
    for(i = 0; i < ncomps; i++)
    {
        char tmp[20], *kc = NULL;
        FmsComponent comp;
        if(FmsMeshGetComponent(mesh, i, &comp))
            E_RETURN(7);

        if(!comp)
            E_RETURN(8);

        sprintf(tmp, "%d", (int)i);
        kc = join_keys(kcomps, tmp);
        FmsIOWriteComponent(ctx, io, kc, comp);
        FREE(kc);
    }
    FREE(kcomps);

    char *ktags = join_keys(key, "Tags");
    for(i = 0; i < ntags; i++)
    {
        char tmp[20], *kt;
        FmsTag t;
        if(FmsMeshGetTag(mesh, i, &t))
            E_RETURN(9);

        sprintf(tmp, "%d", (int)i);
        kt = join_keys(ktags, tmp);
        FmsIOWriteTag(ctx, io, kt, t);
        FREE(kt);
    }
    FREE(ktags);

    return 0;
}

static int
FmsIOReadFmsMesh(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsIOMeshInfo *mesh_info) {\
    int err = 0;
    FmsInt i = 0;
    if(!mesh_info) E_RETURN(1);

    // Get partition info
    {
        char *kpinfo = join_keys(key, "PartitionInfo");
        FmsIntType type;
        FmsInt n = 0;
        err = (*io->get_typed_int_array)(ctx, kpinfo, &type, (void*)&mesh_info->partition_info, &n);
        FREE(kpinfo);
        if(err)
            E_RETURN(2);
    }

    // Get num domain names
    {
        char *kndnames = join_keys(key, "NumDomainNames");
        err = (*io->get_int)(ctx, kndnames, &mesh_info->ndomain_names);
        FREE(kndnames);
        if(err)
            E_RETURN(3);
    }

    // Get num components
    {
        char *kncomps = join_keys(key, "NumComponents");
        err = (*io->get_int)(ctx, kncomps, &mesh_info->ncomponents);
        FREE(kncomps);
        if(err)
            E_RETURN(4);
    }

     // Get num tags
     {
         char *kntags = join_keys(key, "NumTags");
         err = (*io->get_int)(ctx, kntags, &mesh_info->ntags);
         FREE(kntags);
         if(err)
            E_RETURN(5);
     }

    // Get domain names
    if(mesh_info->ndomain_names)
    {
        char *kdom_names = join_keys(key, "DomainNames");
        mesh_info->domain_names = calloc(sizeof(FmsIODomainNameInfo), mesh_info->ndomain_names);
        for(i = 0; i < mesh_info->ndomain_names; i++)
        {
            char temp[21], *dk = NULL;
            sprintf(temp, FMS_LU, i);
            dk = join_keys(kdom_names, temp);
            err |= FmsIOReadFmsDomainName(ctx, io, dk, &mesh_info->domain_names[i]);
            FREE(dk);
        }
        FREE(kdom_names);
        if(err)
            E_RETURN(6);
    }

    // Get components
    if(mesh_info->ncomponents)
    {
        char *kcomps = join_keys(key, "Components");
        mesh_info->components = calloc(sizeof(FmsIOComponentInfo), mesh_info->ncomponents);
        for(i = 0; i < mesh_info->ncomponents; i++) 
        {
            char temp[21], *kc = NULL;
            sprintf(temp, FMS_LU, i);
            kc = join_keys(kcomps, temp);
            err |= FmsIOReadFmsComponent(ctx, io, kc, &mesh_info->components[i]);
            FREE(kc);
        }
        FREE(kcomps);
        if(err)
            E_RETURN(7);
    }

    // Get tags
    if(mesh_info->ntags)
    {
        char *ktags = join_keys(key, "Tags");
        mesh_info->tags = calloc(sizeof(FmsIOTagInfo), mesh_info->ntags);
        for(i = 0; i < mesh_info->ntags; i++)
        {
            char temp[21], *kt = NULL;
            sprintf(temp, FMS_LU, i);
            kt = join_keys(ktags, temp);
            err |= FmsIOReadFmsTag(ctx, io, kt, &mesh_info->tags[i]);
            FREE(kt);
        }
        FREE(ktags);
        if(err)
            E_RETURN(8);   
    }
    return 0;
}

/**
@brief Write FmsFieldDescriptor to the output I/O context.
@param ctx The context
@param io  The I/O functions that operate on the context.
@param key The key that identifies the data in the output.
@param fd  The FMS field descriptor.
@return 0 on success, non-zero otherwise.
*/
static int
FmsIOWriteFmsFieldDescriptor(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, 
    FmsFieldDescriptor fd)
{
    // Write the fd's name
    const char *fd_name = NULL;
    if(FmsFieldDescriptorGetName(fd, &fd_name))
        E_RETURN(1);

    char *key_name = join_keys(key, "Name");
    (*io->add_string)(ctx, key_name, fd_name);
    FREE(key_name);

    // Write which compoenent the fd refers to    
    FmsComponent fc = NULL;
    if(FmsFieldDescriptorGetComponent(fd, &fc))
        E_RETURN(2);

    // Q: Do FieldDescriptors ALWAYS refer to a component?
    if(!fc)
        E_RETURN(3);

    const char *fc_name = NULL;
    if(FmsComponentGetName(fc, &fc_name))
        E_RETURN(4);

    char *key_fc_name = join_keys(key, "ComponentName");
    (*io->add_string)(ctx, key_fc_name, fc_name);
    FREE(key_fc_name);
    
    // Write fd type
    FmsFieldDescriptorType fd_type;
    if(FmsFieldDescriptorGetType(fd, &fd_type))
    {
        E_RETURN(5);
    }

    char *kfd_type = join_keys(key, "Type");
    (*io->add_int)(ctx, kfd_type, (int)fd_type);
    FREE(kfd_type);

    // Write fd fixed order
    FmsFieldType field_type;
    FmsBasisType basis_type;
    FmsInt fo[3];
    if(FmsFieldDescriptorGetFixedOrder(fd, &field_type, &basis_type, &fo[2]))
        E_RETURN(6);
    fo[0] = (FmsInt)field_type; fo[1] = (FmsInt)basis_type;

    char *kfield_type = join_keys(key, "FixedOrder");
    (*io->add_int_array)(ctx, kfield_type, fo, 3);
    FREE(kfield_type);

    // Write fd num dofs
    FmsInt ndofs = 0;
    if(FmsFieldDescriptorGetNumDofs(fd, &ndofs))
        E_RETURN(7);
    
    char *kndofs = join_keys(key, "NumDofs");
    (*io->add_int)(ctx, kndofs, ndofs);
    FREE(kndofs);
    return 0;
}

static int
FmsIOReadFmsFieldDescriptor(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsIOFieldDescriptorInfo *fd_info)
{
    int err = 0;

    // Get the fd's name
    {
        char *kfd_name = join_keys(key, "Name");
        err = (*io->get_string)(ctx, kfd_name, &fd_info->name);
        FREE(kfd_name);
        if(err)
            E_RETURN(1);
    }

    // Get which compoenent the fd refers to   
    {
        char *kcomp_name = join_keys(key, "ComponentName");
        err = (*io->get_string)(ctx, kcomp_name, &fd_info->component_name);
        FREE(kcomp_name);
        if(err)
            E_RETURN(2);
    }

    // Get the fd type
    {
        char *kfd_type = join_keys(key, "Type");
        FmsInt fdt = 0;
        err = (*io->get_int)(ctx, kfd_type, &fdt);
        FREE(kfd_type);
        if(err)
            E_RETURN(3);
        fd_info->type = (FmsFieldDescriptorType)fdt;
    }
    // (There is not setter for fd_type, might not have to write or read)

    // Get FixedOrder
    {
        char *kfo = join_keys(key, "FixedOrder");
        FmsIntType type;
        FmsInt n = 0;
        err = (*io->get_typed_int_array)(ctx, kfo, &type, (void*)fd_info->fixed_order, &n);
        FREE(kfo);
        // We know is is an array of FmsInts that is length 3
        if(err || type != FMS_UINT64 || n != 3u)
            E_RETURN(4);
    }

    // Get ndofs, no setter either, might not have to write or read
    {
        char *kndofs = join_keys(key, "NumDofs");
        err = (*io->get_int)(ctx, kndofs, &fd_info->num_dofs);
        FREE(kndofs);
        if(err)
            E_RETURN(5);
    }

    return 0;
}

/**
@brief Write FmsField to the output I/O context.
@param ctx The context
@param io  The I/O functions that operate on the context.
@param key The key that identifies the data in the output.
@param field  The FMS field.
@return 0 on success, non-zero otherwise.
*/
static int
FmsIOWriteFmsField(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsField field)
{
    const char *field_name;
    if(FmsFieldGetName(field, &field_name))
        E_RETURN(1);
    
    char *kfield_name = join_keys(key, "Name");
    (*io->add_string)(ctx, kfield_name, field_name);
    FREE(kfield_name);

    FmsFieldDescriptor fd = NULL;
    FmsInt num_vec_comp =  0;
    FmsLayoutType lt; FmsScalarType dt;
    const void *data;
    if(FmsFieldGet(field, &fd, &num_vec_comp, &lt, &dt, &data))
        E_RETURN(2);
    
    char *klt = join_keys(key, "LayoutType");
    (*io->add_int)(ctx, klt, (FmsInt)lt);
    FREE(klt);

    char *knc = join_keys(key, "NumberOfVectorComponents");
    (*io->add_int)(ctx, knc, num_vec_comp);
    FREE(knc);

    if(!fd)
        E_RETURN(3);

    const char *fdname = NULL;
    if(FmsFieldDescriptorGetName(fd, &fdname))
        E_RETURN(4);
    char *kfdname = join_keys(key, "FieldDescriptorName");
    (*io->add_string)(ctx, kfdname, fdname);
    FREE(kfdname);

    FmsMetaData md;
    if(FmsFieldGetMetaData(field, &md))
        E_RETURN(4);
    
    if(md)
    {
        char *kmd = join_keys(key, "MetaData");
        FmsIOWriteFmsMetaData(ctx, io, kmd, md);
        FREE(kmd);
    }

    FmsInt ndofs = 0;
    FmsFieldDescriptorGetNumDofs(fd, &ndofs);

    char *kdata = join_keys(key, "Data");
    (*io->add_scalar_array)(ctx, kdata, dt, data, num_vec_comp * ndofs);
    FREE(kdata);


    /* NOTES:
        1. What is in the MetaData?
        2. If each field points to a field descriptor - do should I write FieldDescriptors here?
            * Maybe multiple fields can point to the same FieldDescriptor in which case 
                I should just write them seperate and refer to them. */

    return 0;
}

static int
FmsIOReadFmsField(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsIOFieldInfo *field_info) {
    int err = 0;

    // Get Field name
    {
        char *kfield_name = join_keys(key, "Name");
        err = (*io->get_string)(ctx, kfield_name, &field_info->name);
        FREE(kfield_name);
        if(err)
            E_RETURN(1);
    }
    
    // Get layout type
    {
        char *klt = join_keys(key, "LayoutType");
        FmsInt lt = 0;
        err = (*io->get_int)(ctx, klt, &lt);
        FREE(klt);
        if(err)
            E_RETURN(2);
        field_info->layout = (FmsLayoutType)lt;
    }

    // Get Num vec components
    {
        char *knum_vec_comps = join_keys(key, "NumberOfVectorComponents");
        err = (*io->get_int)(ctx, knum_vec_comps, &field_info->num_vec_comps);
        FREE(knum_vec_comps);
        if(err)
            E_RETURN(3);
    }

    // Get field descriptor name
    {
        char *kfd_name = join_keys(key, "FieldDescriptorName");
        err = (*io->get_string)(ctx, kfd_name, &field_info->fd_name);
        FREE(kfd_name);
        if(err)
            E_RETURN(4);
    }

    // Get scalar data
    {
        char *kdata = join_keys(key, "Data");
        err = (*io->get_scalar_array)(ctx, kdata, &field_info->data_type, &field_info->data, &field_info->data_size);
        FREE(kdata);
        if(err)
            E_RETURN(5);
    }
    
    // Check for MetaData
    {
        char *kmd = join_keys(key, "MetaData");
        if((*io->has_path)(ctx, kmd))
        {
            /* TODO: Put the function call here when it's done
            err = FmsIOReadFmsMetaData(ctx, io, kmd, md); */
            err = 1;
        }
        FREE(kmd);
        if(err)
            E_RETURN(6); 
    }
    return 0;
}

/**
@brief Write FmsDataCollection to the output I/O context.
@param ctx The context
@param io  The I/O functions that operate on the context.
@param key The key that identifies the data in the output.
@param dc  The FMS data collection.
@return 0 on success, non-zero otherwise.
*/
static int
FmsIOWriteFmsDataCollection(FmsIOContext *ctx, FmsIOFunctions *io, const char *key,
    FmsDataCollection dc)
{
    int err = 0;
    FmsMesh mesh;
    FmsFieldDescriptor *fds = NULL;
    FmsField *fields;
    FmsMetaData md;
    FmsInt num_fds = 0, num_fields = 0;
    FmsInt i;
    char *name = NULL;

    /** NOTE: FmsDataCollection does not provide a method to get the name. We need one. */
    if(FmsDataCollectionGetName(dc, &name) == 0)
    {
        char *n_key = join_keys(key, "Name");
        err = (*io->add_string)(ctx, n_key, name);
        if(err)
            E_RETURN(1);
    }

    if(FmsDataCollectionGetFieldDescriptors(dc, &fds, &num_fds) == 0)
    {
        char *fd_key = NULL, *nfds_key = NULL;

        nfds_key = join_keys(key, "NumberOfFieldDescriptors");
        if((*io->add_int)(ctx, nfds_key, num_fds) != 0)
        {
            FREE(nfds_key);
            E_RETURN(2);
        }
        FREE(nfds_key);

        fd_key = join_keys(key, "FieldDescriptors");
        for(i = 0; i < num_fds && err == 0; ++i)
        {
             char tmp[20], *fdi_key = NULL;
             sprintf(tmp, "%d", (int)i);
             fdi_key = join_keys(fd_key, tmp);

             /* Writing DataCollection/FieldDescriptors/0/... */
             err = FmsIOWriteFmsFieldDescriptor(ctx, io, fdi_key, fds[i]);

             FREE(fdi_key);
        }

        FREE(fd_key);
        if(err)
            E_RETURN(3);
    }
    else
    {
        E_RETURN(4);
    }

    if(FmsDataCollectionGetFields(dc, &fields, &num_fields) == 0)
    {
        char *f_key = NULL, *nfs_key = NULL;
        nfs_key = join_keys(key, "NumberOfFields");
        if((*io->add_int)(ctx, nfs_key, num_fields) != 0)
        {
            FREE(nfs_key);
            E_RETURN(5);
        }
        FREE(nfs_key);

        f_key = join_keys(key, "Fields");
        for(i = 0; i < num_fields && err == 0; ++i)
        {
             char tmp[20], *fi_key = NULL;
             sprintf(tmp, "%d", (int)i);
             fi_key = join_keys(f_key, tmp);
             
             err = FmsIOWriteFmsField(ctx, io, fi_key, fields[i]);

             FREE(fi_key);
        }

        FREE(f_key);
        if(err)
            E_RETURN(6);
    }
    else
    {
        E_RETURN(7);
    }

    if(FmsDataCollectionGetMesh(dc, &mesh) == 0)
    {
        char *m_key = join_keys(key, "Mesh");
        err = FmsIOWriteFmsMesh(ctx, io, m_key, mesh);
        FREE(m_key);
        if(err)
            E_RETURN(8);
    }
    else
    {
        E_RETURN(9);
    }

    if(FmsDataCollectionGetMetaData(dc, &md) == 0)
    {
        if(md)
        {
            char *md_key = join_keys(key, "MetaData");
            err = FmsIOWriteFmsMetaData(ctx, io, md_key, md);
            FREE(md_key);
            if(err)
                E_RETURN(10);
        }
    }
    else
    {
        E_RETURN(11);
    }

    return 0;
}

static int
FmsIOReadFmsDataCollection(FmsIOContext *ctx, FmsIOFunctions *io, const char *key,
    FmsIODataCollectionInfo *data_collection)
{
    if(!ctx) E_RETURN(1);
    if(!io) E_RETURN(2);
    if(!key) E_RETURN(3);
    if(!data_collection) E_RETURN(4);
    int err = 0;
    FmsInt i = 0;

    data_collection = calloc(sizeof(FmsIODataCollectionInfo), 1);
    if(!data_collection)
        E_RETURN(5);

    // Process DataCollection's Name
    {
        char *kname = join_keys(key, "Name");
        err = (*io->get_string)(ctx, kname, &data_collection->name);
        FREE(kname);
        if(err)
            E_RETURN(6);
    }

    if(!data_collection->name)
        E_RETURN(7);

    // Process FieldDescriptors
    {
        char *knfds = join_keys(key, "NumberOfFieldDescriptors");
        err = (*io->get_int)(ctx, knfds, &data_collection->nfds);
        FREE(knfds);
        if(err)
            E_RETURN(8);

        // Make sure there were actually any FieldDescriptors
        if(data_collection->nfds)
        {
            char *kfds = join_keys(key, "FieldDescriptors");
            data_collection->fds = calloc(sizeof(FmsIOFieldDescriptorInfo), data_collection->nfds);
            err = 0;
            for(i = 0; i < data_collection->nfds; i++)
            {
                char tmp[20], *kfdi = NULL;
                sprintf(tmp, "%d", (int)i);
                kfdi = join_keys(kfds, tmp);
                err |= FmsIOReadFmsFieldDescriptor(ctx, io, kfdi, &data_collection->fds[i]);
                FREE(kfdi);
            }
            FREE(kfds);
            if(err)
                E_RETURN(11);
        }
    }

    // Process Fields
    {
        char *knfs = join_keys(key, "NumberOfFields");
        err = (*io->get_int)(ctx, knfs, &data_collection->nfields);
        FREE(knfs);
        if(err)
            E_RETURN(12);

        // Make sure there were actually any Fields
        if(data_collection->nfields)
        {
            char *kfs = join_keys(key, "Fields");
            data_collection->fields = calloc(sizeof(FmsIOFieldInfo), data_collection->nfields);
            err = 0;
            for(i = 0; i < data_collection->nfields; i++)
            {
                char tmp[20], *kfi = NULL;
                sprintf(tmp, "%d", (int)i);
                kfi = join_keys(kfs, tmp);
                err |= FmsIOReadFmsField(ctx, io, kfi, &data_collection->fields[i]);
                FREE(kfi);
            }
            FREE(kfs);
            if(err)
                E_RETURN(13);
        }
    }

    // Process mesh
    {
        char *kmesh = join_keys(key, "Mesh");
        err = FmsIOReadFmsMesh(ctx, io, kmesh, &data_collection->mesh_info);
        FREE(kmesh);
        if(err)
            E_RETURN(14);
    }


    // Process MetaData (if we have it)
    {
        char *kmd = join_keys(key, "MetaData");
        if((*io->has_path)(ctx, kmd))
        {
            data_collection->md = calloc(sizeof(FmsIOMetaDataInfo), 1);
            if(data_collection->md)
                err = FmsIOReadFmsMetaData(ctx, io, kmd, data_collection->md);
        }
        FREE(kmd);
        if(err)
            E_RETURN(15);
    }

    return 0;
}

/*****************************************************************************
*** FMS Public Functions
*****************************************************************************/

int
FmsIOWrite(const char *filename, const char *protocol, FmsDataCollection dc)
{
    FmsIOContext ctx;
    FmsIOFunctions io;

    if(filename == NULL) E_RETURN(1);
    if(protocol == NULL) protocol = "ascii";

#ifdef FMS_HAVE_CONDUIT
    /* Conduit has a function to enumerate its protocols that it supports. Use that later. */
    if(strcmp(protocol, "json") == 0 ||
       strcmp(protocol, "yaml") == 0 ||
       strcmp(protocol, "hdf5") == 0 ||
       strcmp(protocol, "silo") == 0 ||
       strcmp(protocol, "adios") == 0)
    {
        FmsIOContextInitializeConduit(&ctx, protocol);
        FmsIOFunctionsInitializeConduit(&io);
    }
    else
    {
        FmsIOContextInitialize(&ctx);
        FmsIOFunctionsInitialize(&io);
    }
#else
    FmsIOContextInitialize(&ctx);
    FmsIOFunctionsInitialize(&io);
#endif

    /* "open" the file. */
    if((*io.open)(&ctx, filename, "wt") == 0)
    {
        if((*io.add_int)(&ctx, "FMS", (int)FMS_INTERFACE_VERSION) != 0)
        {
            /* "close" the file. */
            (*io.close)(&ctx);
            E_RETURN(4);
        }

        if(FmsIOWriteFmsDataCollection(&ctx, &io, "DataCollection", dc) != 0)
        {
            /* "close" the file. */
            (*io.close)(&ctx);
            E_RETURN(5);
        }

        if((*io.close)(&ctx) != 0)
            E_RETURN(6);
    }
    else
    {
        E_RETURN(3);
    }

    return 0;
}

int
FmsIORead(const char *filename, const char *protocol, FmsDataCollection *dc)
{
    FmsIOContext ctx;
    FmsIOFunctions io;

    if(!filename) E_RETURN(1);
    if(!protocol) protocol = "ascii";
    if(!dc) E_RETURN(3);

#ifdef FMS_HAVE_CONDUIT
    /* Conduit has a function to enumerate its protocols that it supports. Use that later. */
    if(strcmp(protocol, "json") == 0 ||
       strcmp(protocol, "yaml") == 0 ||
       strcmp(protocol, "hdf5") == 0 ||
       strcmp(protocol, "silo") == 0 ||
       strcmp(protocol, "adios") == 0)
    {
        FmsIOContextInitializeConduit(&ctx, protocol);
        FmsIOFunctionsInitializeConduit(&io);
    }
    else
    {
        FmsIOContextInitialize(&ctx);
        FmsIOFunctionsInitialize(&io);
    }
#else
    FmsIOContextInitialize(&ctx);
    FmsIOFunctionsInitialize(&io);
#endif

    /* "open" the file. */
    if((*io.open)(&ctx, filename, "r") == 0)
    {
        int err = 0;
        FmsInt version = 0;
        FmsIODataCollectionInfo dc_info;
        err = (*io.get_int)(&ctx, "FMS", &version);
        if(err != 0 || version != (int)FMS_INTERFACE_VERSION)
        {
            /* "close" the file. */
            (*io.close)(&ctx);
            E_RETURN(5);
        }

        if(FmsIOReadFmsDataCollection(&ctx, &io, "DataCollection", &dc_info) != 0)
        {
            /* "close" the file. */
            /* TODO: Destroy the FmsIODataCollectionInfo */
            (*io.close)(&ctx);
            E_RETURN(6);
        }

        if((*io.close)(&ctx) != 0)
            E_RETURN(7);
    }
    else
    {
        E_RETURN(4);
    }

    return 0;
}

