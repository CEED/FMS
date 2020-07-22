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
#include <fmsio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONDUIT
#include <conduit.h>
#endif

#define FREE(PTR) if((PTR) != NULL) free(PTR)

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
#ifdef HAVE_CONDUIT
    /* Conduit C-API*/
    conduit_node *root;

#ifdef MEANDERING_THOUGHTS
    conduit_node *stack[100];
    int           stack_top;
#endif
    char *filename;
    char *protocol;
    int   writing;
    
#else
    FILE *fp;
#endif
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
    int (*add_int_array)(FmsIOContext *ctx, const char *path, const FmsInt *values, size_t n);
    int (*add_float)(FmsIOContext *ctx, const char *path, float value);
    int (*add_float_array)(FmsIOContext *ctx, const char *path, const float *values);
    int (*add_double)(FmsIOContext *ctx, const char *path, double value);
    int (*add_double_array)(FmsIOContext *ctx, const char *path, const double *value, size_t n);
    int (*add_string)(FmsIOContext *ctx, const char *path, const char *value);
    int (*add_string_array)(FmsIOContext *ctx, const char *path, const char **value, size_t n);

    /* TODO: functions for reading path+value*/
    int (*get_int)(FmsIOContext *ctx, const char *path, int *value);
    /*...*/

} FmsIOFunctions;

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

static int
FmsIOAddInt(FmsIOContext *ctx, const char *path, FmsInt value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    fprintf(ctx->fp, "%s: %llu\n", path, value);
    return 0;
}

static int
FmsIOAddIntArray(FmsIOContext *ctx, const char *path, const FmsInt *values, size_t n)
{
    size_t i;
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
#if 1
    /* Should we make it YAML-like?*/
    fprintf(ctx->fp, "%s: [", path);
    for(i = 0; i < n; ++i)
    {
        fprintf(ctx->fp, "%llu", values[i]);
        if(i < n-1)
            fprintf(ctx->fp, ", ");
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
    fprintf(ctx->fp, "%s: %lg\n", path, value);
    return 0;
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

int
FmsIOGetInt(FmsIOContext *ctx, const char *path, int *value)
{
    char *line = NULL;
    ssize_t nread = 0;
    size_t len = 0;
    int retval = 0;
    if(ctx) E_RETURN(1);
    if(path) E_RETURN(2);
    if(value) E_RETURN(3);
    if((nread = getline(&line, &len, ctx->fp)) != -1)
    {
        size_t klen = strlen(path);
        if(strncmp(line, path, klen) == 0)
        {
            if(sscanf(line+klen+1, "%d", value) != 1)
                retval = 4;
        }
        else
            retval = 5;
        FREE(line);
    }
    return retval;
}

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
        obj->add_float = FmsIOAddFloat;
        obj->add_double = FmsIOAddDouble;
        obj->add_string = FmsIOAddString;
        /*...*/
        obj->get_int = FmsIOGetInt;
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
#ifdef HAVE_CONDUIT
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
FMSIOCloseConduit(FmsIOContext *ctx)
{
    if(!fp) E_RETURN(1);

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
FmsIOAddIntConduit(FmsIOContext *ctx, const char *path, FmsInt value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    /*conduit_node_set_path_int(ctx->root, path, value);*/
    conduit_node_set_path_unsigned_long(ctx->root, path, value);
    return 0;
}

static int
FmsIOAddIntArrayConduit(FmsIOContext *ctx, const char *path, const int *values, size_t n)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    if(!values) E_RETURN(3);
    /*conduit_node_set_path_int_ptr(ctx->root, path, (int *)values, n);*/
    conduit_node_set_path_uint64_ptr(ctx->root, path, values, n);
    return 0;
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
FmsIOAddStringConduit(FmsIOContext *ctx, const char *path, const char *value)
{
    if(!ctx) E_RETURN(1);
    if(!path) E_RETURN(2);
    /*conduit_node_set_path_char8_str(ctx->root, path, value);*/
    conduit_node_set_path_external_char8_str(ctx->root, path, value);
    return 0;
}

int
FmsIOGetIntConduit(FmsIOContext *ctx, const char *path, int *value)
{
    conduit_node *node = NULL;
    if(ctx) E_RETURN(1);
    if(path) E_RETURN(2);
    if(value) E_RETURN(3);
    if((node = conduit_node_fetch(ctx->root, path)) != NULL)
    {
        const conduit_datatype *dt = conduit_node_dtype(node);
        if(!dt)
        {
            if(conduit_datatype_is_number(dt))
            {
                *value = conduit_node_fetch_path_as_int(ctx->root, path);
            }
            else
            {
                E_RETURN(5);
            }
        }
    }
    else
    {
        E_RETURN(4);
    }
    return 0
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
        for(i = 0; i < 100; ++i)
            ctx->stack[i] = NULL;
        ctx->stack_top = -1;

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
        obj->add_float = FmsIOAddFloatConduit;
        obj->add_double = FmsIOAddDoubleConduit;
        obj->add_string = FmsIOAddStringConduit;
        /*...*/
        obj->get_int = FmsIOGetIntConduit;
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
    /** TODO: write me */
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
    /** TODO: write me */
#if 0
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
#endif

    
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
    /** TODO: write me */
    return 0;
}

/**
@brief Write FmsMetaData to the output I/O context.
@param ctx The context
@param io  The I/O functions that operate on the context.
@param key The key that identifies the data in the output.
@param md  The FMS metadata.
@return 0 on success, non-zero otherwise.
*/
static int
FmsIOWriteFmsMetaData(FmsIOContext *ctx, FmsIOFunctions *io, const char *key, FmsMetaData md)
{
    /** TODO: write me */
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

    if(FmsDataCollectionGetMesh(dc, &mesh) == 0)
    {
        char *m_key = join_keys(key, "Mesh");
        err = FmsIOWriteFmsMesh(ctx, io, m_key, mesh);
        FREE(m_key);
        if(err)
            E_RETURN(1);
    }
    else
    {
        E_RETURN(2);
    }

    if(FmsDataCollectionGetFieldDescriptors(dc, &fds, &num_fds) == 0)
    {
        char *fd_key = NULL, *fdlen_key = NULL;
        fd_key = join_keys(key, "FieldDescriptors");

        fdlen_key = join_keys(fd_key, "Length");
        if((*io->add_int)(ctx, fd_key, num_fds) != 0)
        {
            FREE(fd_key);
            FREE(fdlen_key);
            E_RETURN(3);
        }
        FREE(fdlen_key);

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
            E_RETURN(4);
    }
    else
    {
        E_RETURN(5);
    }

    if(FmsDataCollectionGetFields(dc, &fields, &num_fields) == 0)
    {
        char *f_key = NULL, *flen_key = NULL;
        f_key = join_keys(key, "Fields");

        flen_key = join_keys(f_key, "Length");
        if((*io->add_int)(ctx, f_key, num_fields) != 0)
        {
            FREE(f_key);
            FREE(flen_key);
            E_RETURN(6);
        }
        FREE(flen_key);

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
            E_RETURN(7);
    }
    else
    {
        E_RETURN(8);
    }

    if(FmsDataCollectionGetMetaData(dc, &md) == 0)
    {
        char *md_key = join_keys(key, "MetaData");
        err = FmsIOWriteFmsMetaData(ctx, io, md_key, md);
        FREE(md_key);
        if(err)
            E_RETURN(9);
    }
    else
    {
        E_RETURN(10);
    }

    return 0;
}

static int
FmsIOReadFmsDataCollection(FmsIOContext *ctx, FmsIOFunctions *io, const char *key,
    FmsDataCollection *dc)
{
    /*TODO: write me. We should be able to mostly copy structure of FmsIOWriteFmsDataCollection */
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
    if(protocol == NULL) E_RETURN(2);

#ifdef HAVE_CONDUIT
    /* Conduit has a function to enumerate its protocols that it supports. Use that later. */
    if(strcmp(protocol, "json") == 0 ||
       strcmp(protocol, "yaml") == 0 ||
       strcmp(protocol, "hdf5") == 0 ||
       strcmp(protocol, "silo") == 0 ||
       strcmp(protocol, "adios") == 0)
    {
        FmsIOContextInitializeConduit(&ctx);
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
    if(!protocol) E_RETURN(2);
    if(!dc) E_RETURN(3);

#ifdef HAVE_CONDUIT
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
        int err = 0, version = 0;
        err = (*io.get_int)(&ctx, "FMS", &version);
        if(err != 0 || version != (int)FMS_INTERFACE_VERSION)
        {
            /* "close" the file. */
            (*io.close)(&ctx);
            E_RETURN(5);
        }

        if(FmsIOReadFmsDataCollection(&ctx, &io, "DataCollection", dc) != 0)
        {
            /* "close" the file. */
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

