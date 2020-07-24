// Example of using the FMS interface to describe the following 1-element mesh
//
//  <3>----(2)----<2>
//   |             |
//   |             |
//  (3)    [0]    (1)
//   |             |
//   |             |
//  <0>----(0)----<1>
//
// where <V> = vertex id, (E) = edge id and [Q] = quadrilateral id. The edge and
// element orientation are counter-clockwise:
//
// (0) = <0>,<1> ; (1) = <1>,<2> ; (2) = <2>,<3> ; (3) = <3>,<0>
// [0] = (0), (1), (2), (3)
//
// The mesh has one domain and two components: a volume containing element [0]
// and boundary containing edge (2). The volume component is tagged with the
// "material" tag.

#include <fms.h>

#define TEST_FMSIO 1
#if TEST_FMSIO
#include <fmsio.h>
#endif

int main(int argc, const char *argv[]) {

  // Create the Mesh object
  FmsMesh mesh;
  FmsMeshConstruct(&mesh);
  // FmsMeshSetPartitionId(mesh, 0);

  // Mesh Domains
  FmsDomain *domain;
  FmsMeshAddDomains(mesh, "domains", 1, &domain);

  FmsDomainSetNumVertices(domain[0], 4);
  FmsDomainSetNumEntities(domain[0], FMS_EDGE, FMS_INT32, 4);
  FmsDomainSetNumEntities(domain[0], FMS_QUADRILATERAL, FMS_INT32, 1);

  // Edges
  const int edge_vert[] = {0,1, 1,2, 2,3, 3,0};
  FmsDomainAddEntities(domain[0], FMS_EDGE, NULL, FMS_INT32, edge_vert, 4);

  // Faces
  const int quad_edge[] = {0, 1, 2, 3};
  FmsDomainAddEntities(domain[0], FMS_QUADRILATERAL, NULL,
                       FMS_INT32, quad_edge, 1);

  // Mesh Components
  FmsComponent volume;
  FmsMeshAddComponent(mesh, "volume", &volume);
  FmsComponentAddDomain(volume, domain[0]);

  FmsComponent boundary;
  FmsMeshAddComponent(mesh, "boundary", &boundary);
  FmsInt part_id;
  FmsComponentAddPart(boundary, domain[0], &part_id);
  const int bdr_edge[] = {2};
  FmsComponentAddPartEntities(boundary, part_id, FMS_EDGE, FMS_INT32,
                              FMS_INT32, FMS_INT32, NULL, bdr_edge, NULL, 1);

  FmsComponentAddRelation(volume, 1); // 1 = index of "boundary" component

  // Mesh Tags
  FmsTag material;
  FmsMeshAddTag(mesh, "material", &material);
  FmsTagSetComponent(material, volume);
  const int material_tags[] = {1};
  FmsTagSet(material, FMS_INT32, FMS_INT32, material_tags, 1);

  int tvs[1] = {1};
  const char * const descriptions[] = {"Test"};
  FmsTagAddDescriptions(material, FMS_INT32, tvs, descriptions, 1);

  // Finalize the construction of the Mesh object
  FmsMeshFinalize(mesh);
  // Perform some consistency checks.
  FmsMeshValidate(mesh);

  // Coordinates Field
  // Define data collection
  FmsDataCollection dc;
  FmsDataCollectionCreate(mesh, "data collection", &dc);

  // Define field descriptor
  FmsFieldDescriptor coords_fd;
  FmsDataCollectionAddFieldDescriptor(dc, "coords descriptor", &coords_fd);
  FmsFieldDescriptorSetComponent(coords_fd, volume);
  FmsInt coords_order = 1;
  FmsFieldDescriptorSetFixedOrder(coords_fd, FMS_CONTINUOUS,
                                  FMS_NODAL_GAUSS_CLOSED, coords_order);

  // Define the coordinates field
  FmsField coords;
  FmsDataCollectionAddField(dc, "coords", &coords);
  const double coords_data[] = {
    0.,0.,
    1.,0.,
    1.,1.,
    0.,1.
  };
  FmsInt sdim = 2; // A 2D mesh embedded in "sdim"-dimensional space.
  FmsFieldSet(coords, coords_fd, sdim, FMS_BY_VDIM, FMS_DOUBLE, coords_data);

  FmsComponentSetCoordinates(volume, coords);

  // Use the mesh to perform computations
  if(argc > 1)
    FmsIOWrite("simple_mesh.fms", argv[1], dc);
  else
    FmsIOWrite("simple_mesh.fms.ascii", NULL, dc);

  // Destroy the data collection: destroys the FmsMesh and all other linked Fms
  // objects.
  FmsDataCollectionDestroy(&dc);

  return 0;
}