//
// Demo: converting a mesh from PUMI to FMS to MFEM
//

#include <fms.h>
#include <apf.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <apfMDS.h>
#include <apfNumbering.h>
#include <gmi.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <crv.h>
#include <vector>
#include <cassert>
#include <iostream>
#include <limits.h>
#include <mfem.hpp>

#define TEST_FMSIO
#ifdef TEST_FMSIO
#include <fmsio.h>
#endif

void TestMakeFmsQuadMesh(FmsDataCollection *dc_ptr);
void TestMakeFmsTetMesh(FmsDataCollection *dc_ptr);
void TestMakeFmsHexMesh(FmsDataCollection *dc_ptr);

int ConvertFmsToMfem(FmsDataCollection dc, mfem::Mesh **mesh_p);

void TestFmsToMfem(FmsDataCollection dc);


int main(int argc, char *argv[]) {

  // Initialize pumi
  int num_proc, myId;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);
  PCU_Comm_Init();
  gmi_register_mesh();

  // Load pumi mesh
  const char *mesh_file = "../data/pumi/serial/pillbox.smb";
  const char *model_file = "../data/pumi/geom/pillbox.dmg";
  apf::Mesh2* pumi_mesh;
  pumi_mesh = apf::loadMdsMesh(model_file, mesh_file);

  const int geom_order = 2; // 3rd order has issues ...
  if (geom_order != 2) {
    crv::BezierCurver bc(pumi_mesh, geom_order, 2);
    bc.run();
    // crv::InterpolatingCurver ic(pumi_mesh, geom_order);
    // ic.run();
  }

  // The reference face and element entity ordering follows PUMI convention
  // presented in pages 13 & 14 of https://scorec.rpi.edu/pumi/PUMI.pdf.
  //
  //
  //                       <3>--(5)-<2>    [0] Bottom face containing <0><1><2>
  //         <2>            | \   .   \    [1] Front  face containing <0><1><3>
  //        /   \           |   \.     \   [2] Right  face containing <1><2><3>
  //       /     \         (3)  .(4)    \  [3] Left   face containing <0><2><3>
  //     (2)     (1)        | (2)   \   (1)
  //     /   [0]   \        | .       \   \
  //    /           \       |.          \  \
  //  <0>----(0)----<1>    <0>-----(0)----<1>
  //
  // where <V> = vertex id, (E) = edge id and [T] = Tetrahedral id.

  // First number vertices
  apf::Field* apf_field_crd = pumi_mesh->getCoordinateField();
  apf::FieldShape* crd_shape = apf::getShape(apf_field_crd);
  apf::Numbering* vtx_num_local = apf::createNumbering(pumi_mesh,
                                  "VertexNumbering",
                                  crd_shape, 1);

  // Create the Mesh object
  FmsMesh mesh;
  FmsMeshConstruct(&mesh);
  // CeedMeshSetParallelId(mesh, 0);

  // Mesh Domains
  FmsDomain *domain;
  FmsInt num_domains = 1;
  FmsMeshAddDomains(mesh, "domains", num_domains, &domain);

  // Count and set entities
  FmsInt num_verts = 0;
  apf::MeshIterator* itr = pumi_mesh->begin(0);
  apf::MeshEntity* ent;
  while ((ent = pumi_mesh->iterate(itr))) {
    // IDs start from 0
    apf::number(vtx_num_local, ent, 0, 0, num_verts);
    num_verts++;
  }
  pumi_mesh->end(itr);
  // 0d-entities
  FmsDomainSetNumVertices(domain[0], num_verts);

  FmsInt dim = pumi_mesh->getDimension();
  for (FmsInt n_dim = 1; n_dim <= dim; n_dim++) {
    FmsInt num_ents = countOwned(pumi_mesh, n_dim);
    switch(n_dim) {
    case 1: // 1d-entities
      FmsDomainSetNumEntities(domain[0], FMS_EDGE, FMS_UINT32, num_ents);
      std :: cout << " num FMS_EDGE : " << num_ents << std :: endl;
      break;
    case 2: // 2d-entities
      FmsDomainSetNumEntities(domain[0], FMS_TRIANGLE, FMS_UINT32, num_ents);
      std :: cout << " num FMS_TRIANGLE : " << num_ents << std :: endl;
      break;
    case 3: // 3d-entities
      FmsDomainSetNumEntities(domain[0], FMS_TETRAHEDRON, FMS_UINT32, num_ents);
      std :: cout << " num FMS_TETRAHEDRON : " << num_ents << std :: endl;
      break;
    }
  }

  // Edges
  FmsInt nvtxs = 2;
  FmsInt num_edges = 0;
  int edge[nvtxs];
  itr = pumi_mesh->begin(1);
  while((ent = pumi_mesh->iterate(itr))) {
    // Get vertices
    apf::Downward verts;
    pumi_mesh->getDownward(ent,0,verts);

    for (FmsInt vt = 0; vt < nvtxs; vt++) {
      edge[vt] = apf::getNumber(vtx_num_local, verts[vt], 0, 0);
    }
    //Create edge
    FmsDomainAddEntities(domain[0], FMS_EDGE, 0, FMS_UINT32, edge, 1);
    num_edges++;
  }
  pumi_mesh->end(itr);

  // Faces
  // First number edges
  apf::FieldShape* edge_shape =apf::getConstant(1);
  apf::Numbering* edge_num_local = apf::createNumbering(pumi_mesh,
                                   "EdgeNumbering",
                                   edge_shape, 1);
  itr = pumi_mesh->begin(1);
  FmsInt edgeId = 0;
  while ((ent = pumi_mesh->iterate(itr))) {
    // Assign as local number
    apf::number(edge_num_local, ent, 0, 0, edgeId);
    edgeId++;
  }
  pumi_mesh->end(itr);

  // Define Face based on edge Id's
  FmsInt nedges = 3;
  int face[nedges];
  int num_triangles = 0;
  itr = pumi_mesh->begin(2);
  while ((ent = pumi_mesh->iterate(itr))) {
    // Get edges
    apf::Downward edges;
    pumi_mesh->getDownward(ent, apf::Mesh::EDGE, edges);

    for (FmsInt ed = 0; ed < nedges; ed++) {
      face[ed] = apf::getNumber(edge_num_local, edges[ed], 0, 0);
    }
    // Create face
    FmsDomainAddEntities(domain[0], FMS_TRIANGLE, 0, FMS_UINT32, face, 1);
    num_triangles++;
  }
  pumi_mesh->end(itr);

  // Elements
  // First number face
  apf::FieldShape* face_shape =apf::getConstant(2);
  apf::Numbering* face_num_local = apf::createNumbering(pumi_mesh,
                                   "FaceNumbering",
                                   face_shape, 1);
  itr = pumi_mesh->begin(2);
  FmsInt faceId = 0;
  while ((ent = pumi_mesh->iterate(itr))) {
    // Assign as local number
    apf::number(face_num_local, ent, 0, 0, faceId);
    faceId++;
  }
  pumi_mesh->end(itr);

  // Define Elem based on face Id's
  FmsInt nfaces = 4;
  FmsEntityReordering EntReord;
  int tetOrd[4] = {0,1,3,2}; // PUMI order for tet element faces
  EntReord[FMS_TETRAHEDRON] = &tetOrd[0];
  int elem[nfaces];
  FmsInt num_elems = 0;
  itr = pumi_mesh->begin(3);
  while ((ent = pumi_mesh->iterate(itr))) {
    // Get faces
    apf::Downward faces;
    pumi_mesh->getDownward(ent, apf::Mesh::TRIANGLE, faces);

    for (FmsInt fc = 0; fc < nfaces; fc++) {
      elem[fc] = apf::getNumber(face_num_local, faces[fc], 0, 0);
    }
    // Create elem
    FmsDomainAddEntities(domain[0], FMS_TETRAHEDRON, EntReord,
                         FMS_UINT32, elem, 1);
    num_elems++;
  }
  pumi_mesh->end(itr);


  // Mesh Components
  FmsComponent volume;
  FmsMeshAddComponent(mesh, "volume", &volume);
  FmsComponentAddDomain(volume, domain[0]);

  // Boundary component
  FmsComponent boundary;
  FmsMeshAddComponent(mesh, "boundary", &boundary);
  std::vector<FmsInt> bdr_ents;
  std::vector<FmsInt> bdr_tags;

  itr = pumi_mesh->begin(dim-1);
  FmsInt bdr_cnt = 0;
  while ((ent = pumi_mesh->iterate(itr))) {
    apf::ModelEntity *me = pumi_mesh->toModel(ent);
    if (pumi_mesh->getModelType(me) == (int)(dim-1)) {
      int tag = pumi_mesh->getModelTag(me);
      bdr_ents.push_back(bdr_cnt);
      bdr_tags.push_back(tag);
      bdr_cnt++;
    }
  }
  pumi_mesh->end(itr);

  FmsInt part = 0;
  FmsComponentAddPart(boundary, domain[0], &part);
  FmsComponentAddPartEntities(boundary, part, FMS_TRIANGLE, FMS_UINT32,
                              FMS_INT_TYPE, FMS_INT_TYPE,
                              NULL, &bdr_ents[0], NULL, bdr_ents.size());
  // Add relation 'volume' -> 'boundary'
  FmsComponentAddRelation(volume, 1); // 1 = index of 'boundary' component

  // Boundary tag
  FmsTag tag;
  FmsMeshAddTag(mesh, "boundary", &tag);
  FmsTagSetComponent(tag, boundary);
  FmsTagSet(tag, FMS_INT32, FMS_INT_TYPE, &bdr_tags[0], bdr_tags.size());

  // Add coordinates field
  FmsDataCollection CrdFieldDataCol;
  FmsDataCollectionCreate(mesh, "DataCrd", &CrdFieldDataCol);

  FmsFieldDescriptor CrdFieldDesc;
  FmsDataCollectionAddFieldDescriptor(CrdFieldDataCol, "FieldDesc",
                                      &CrdFieldDesc);
  FmsFieldDescriptorSetComponent(CrdFieldDesc, volume);
  FmsFieldDescriptorSetFixedOrder(CrdFieldDesc, FMS_CONTINUOUS,
                                  FMS_NODAL_GAUSS_CLOSED, geom_order);

  FmsField FmsCrdField;
  FmsInt numVecComp = 3; // Always 3D coordinates

  // Count 0d entities classified on all entities (vertices + nodes)
  apf::Field* pumiCrdField = pumi_mesh->getCoordinateField();
  apf::FieldShape* CrdFieldShape = getShape(pumiCrdField);
  FmsInt numZeroDEntities = num_verts;
  if (CrdFieldShape->hasNodesIn(apf::Mesh::EDGE)) {
    const int nn = CrdFieldShape->countNodesOn(apf::Mesh::EDGE);
    numZeroDEntities += num_edges * FmsInt(nn);
  }
  if (CrdFieldShape->hasNodesIn(apf::Mesh::TRIANGLE)) {
    const int nn = CrdFieldShape->countNodesOn(apf::Mesh::TRIANGLE);
    numZeroDEntities += num_triangles * FmsInt(nn);
  }
  std::cout << "numZeroDEntities : " << numZeroDEntities << std::endl;
  std::vector<double> allPumiCrds(numZeroDEntities * numVecComp);

  // Add vertices coords
  itr = pumi_mesh->begin(0);
  FmsInt counter = 0;
  while ((ent = pumi_mesh->iterate(itr))) {
    unsigned int id = apf::getNumber(vtx_num_local, ent, 0, 0);
    apf::Vector3 Crds;
    pumi_mesh->getPoint(ent,0,Crds);
    int index = id * numVecComp;
    for (int ii=0; ii<(int)numVecComp; ii++) {
      allPumiCrds[index + ii] = Crds[ii];
    }
    counter++;
  }
  pumi_mesh->end(itr);

  // Nodes-classified-on-edge coordinates
  if (CrdFieldShape->hasNodesIn(apf::Mesh::EDGE)) {
    itr = pumi_mesh->begin(1);
    while ((ent = pumi_mesh->iterate(itr))) {
      apf::Vector3 Crds;
      for (int ni = 0; ni < geom_order-1; ni++) {
        pumi_mesh->getPoint(ent,ni,Crds);
        int index = counter * numVecComp;
        for (int ii=0; ii<(int)numVecComp; ii++) {
          allPumiCrds[index + ii] = Crds[ii];
        }
        counter++;
      }
    }
    pumi_mesh->end(itr);
  }

  // Nodes-classified-on-triangle coordinates
  int tri_nodes = ((geom_order-1)*(geom_order-2))/2;
  if (CrdFieldShape->hasNodesIn(apf::Mesh::TRIANGLE)) {
    itr = pumi_mesh->begin(2);
    while ((ent = pumi_mesh->iterate(itr))) {
      apf::Vector3 Crds;
      for (int ni = 0; ni < tri_nodes; ni++) {
        pumi_mesh->getPoint(ent,ni,Crds);
        int index = counter * numVecComp;
        for (int ii=0; ii<(int)numVecComp; ii++) {
          allPumiCrds[index + ii] = Crds[ii];
        }
        counter++;
      }
    }
    pumi_mesh->end(itr);
  }

  // Check all coordinates are added
  assert(counter == numZeroDEntities);

  // Add a coordinate field to the data collection (initializes FmsCrdField)
  FmsDataCollectionAddField(CrdFieldDataCol, "CrdField", &FmsCrdField);

  // Set FMS coordinates field data
  FmsFieldSet(FmsCrdField, CrdFieldDesc, numVecComp, FMS_BY_VDIM, FMS_DOUBLE,
              &allPumiCrds[0]);

  // Set 'FmsCrdField' as the coordinates field for the 'volume' component.
  FmsComponentSetCoordinates(volume, FmsCrdField);


  // Mesh Tags
  FmsTag material;
  FmsMeshAddTag(mesh, "material", &material);
  FmsTagSetComponent(material, volume);
  std::vector<int> ent_tags(num_elems, 314);
  FmsTagSet(material, FMS_UINT32, FMS_UINT32, &ent_tags[0], num_elems);

  // Finalize the construction of the Mesh object
  FmsMeshFinalize(mesh);

  // Perform some consistency checks.
  FmsMeshValidate(mesh);

  FmsDataCollection quad_dc, tet_dc, hex_dc;
  TestMakeFmsQuadMesh(&quad_dc);
  TestMakeFmsTetMesh(&tet_dc);
  TestMakeFmsHexMesh(&hex_dc);

  // Use the data collection
  // TestFmsToMfem(quad_dc);
  // TestFmsToMfem(tet_dc);
  // TestFmsToMfem(hex_dc);
  TestFmsToMfem(CrdFieldDataCol);

#ifdef TEST_FMSIO
  const char *protocol = "ascii";
  FmsIOWrite("quad.fms", protocol, quad_dc);
  FmsIOWrite("tet.fms", protocol, tet_dc);
  FmsIOWrite("hex.fms", protocol, hex_dc);
  FmsIOWrite("CrdFieldDataCol.fms", protocol, CrdFieldDataCol);
#endif

  // Clean-up

  // Destroy the data collection: destroys the mesh and all other linked Fms
  // objects.
  FmsDataCollectionDestroy(&hex_dc);
  FmsDataCollectionDestroy(&tet_dc);
  FmsDataCollectionDestroy(&quad_dc);
  FmsDataCollectionDestroy(&CrdFieldDataCol);

  // Finalize PUMI
  apf::destroyNumbering(vtx_num_local);
  apf::destroyNumbering(edge_num_local);
  apf::destroyNumbering(face_num_local);
  pumi_mesh->destroyNative();
  apf::destroyMesh(pumi_mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}


/* -------------------------------------------------------------------------- */
/* Construct an FMS data collection with one quadrilateral */
/* -------------------------------------------------------------------------- */

void TestMakeFmsQuadMesh(FmsDataCollection *dc_ptr) {
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

  // Finalize the construction of the Mesh object
  FmsMeshFinalize(mesh);
  // Perform some consistency checks.
  FmsMeshValidate(mesh);

  // Coordinates Field
  // Defdine data collection
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

  *dc_ptr = dc;
}


/* -------------------------------------------------------------------------- */
/* Construct an FMS data collection with one tetrahedron */
/* -------------------------------------------------------------------------- */

void TestMakeFmsTetMesh(FmsDataCollection *dc_ptr) {
  FmsMesh mesh;
  FmsMeshConstruct(&mesh);

  FmsDomain *domains;
  FmsMeshAddDomains(mesh, "domain", 1, &domains);
  FmsDomainSetNumVertices(domains[0], 4);

  FmsDomainSetNumEntities(domains[0], FMS_EDGE, FMS_INT32, 6);
  const int edge_vert[] = {0,1, 2,1, 0,2, 0,3, 3,1, 2,3};
  FmsDomainAddEntities(domains[0], FMS_EDGE, NULL, FMS_INT32, edge_vert, 6);

  FmsDomainSetNumEntities(domains[0], FMS_TRIANGLE, FMS_INT32, 4);
  const int face_edge[] = {0,2,1, 0,4,3, 2,3,5, 1,5,4};
  FmsDomainAddEntities(domains[0], FMS_TRIANGLE, NULL, FMS_INT32, face_edge, 4);

  FmsDomainSetNumEntities(domains[0], FMS_TETRAHEDRON, FMS_INT32, 1);
  const int tet_face[] = {0,1,2,3};
  FmsDomainAddEntities(domains[0], FMS_TETRAHEDRON, NULL, FMS_INT32, tet_face,
                       1);

  FmsComponent volume, boundary;
  FmsMeshAddComponent(mesh, "volume", &volume);
  FmsComponentAddDomain(volume, domains[0]);

  FmsMeshAddComponent(mesh, "boundary", &boundary);
  FmsInt part_id;
  FmsComponentAddPart(boundary, domains[0], &part_id);
  const int bdr_face[] = {0,1,2,3};
  FmsComponentAddPartEntities(boundary, part_id, FMS_TRIANGLE, FMS_INT32,
                              FMS_INT32, FMS_INT32, NULL, bdr_face, NULL, 4);

  FmsComponentAddRelation(volume, 1); // 1 = index of "boundary" component

  FmsTag attributes, bdr_attributes;
  FmsMeshAddTag(mesh, "attributes", &attributes);
  FmsTagSetComponent(attributes, volume);
  const int tet_attr[] = {1};
  FmsTagSet(attributes, FMS_UINT8, FMS_INT32, tet_attr, 1);
  const int described_attr[] = {1};
  const char *attr_descr[] = {"material 1"};
  FmsTagAddDescriptions(attributes, FMS_INT32, described_attr, attr_descr, 1);

  FmsMeshAddTag(mesh, "bdr_attributes", &bdr_attributes);
  FmsTagSetComponent(bdr_attributes, boundary);
  const int bdr_attr[] = {1,2,3,4};
  FmsTagSet(bdr_attributes, FMS_UINT8, FMS_INT32, bdr_attr, 4);
  const int bdr_described_attr[] = {1,2,3,4};
  const char *bdr_attr_descr[] = {"face 1", "face 2", "face 3", "face 4"};
  FmsTagAddDescriptions(bdr_attributes, FMS_INT32, bdr_described_attr,
                        bdr_attr_descr, 4);

  FmsMeshFinalize(mesh);
  FmsMeshValidate(mesh);

  FmsDataCollectionCreate(mesh, "data collection", dc_ptr);
  FmsDataCollection dc = *dc_ptr;

  FmsFieldDescriptor coords_fd;
  FmsDataCollectionAddFieldDescriptor(dc, "coords descriptor", &coords_fd);
  FmsFieldDescriptorSetComponent(coords_fd, volume);
  FmsInt order = 1;
  FmsFieldDescriptorSetFixedOrder(coords_fd, FMS_CONTINUOUS,
                                  FMS_NODAL_GAUSS_CLOSED, order);

  FmsField coords;
  FmsDataCollectionAddField(dc, "coords", &coords);
  const double coords_data[] = {
    0.,0.,0.,
    1.,0.,0.,
    0.,1.,0.,
    0.,0.,1.
  };
  FmsFieldSet(coords, coords_fd, 3, FMS_BY_VDIM, FMS_DOUBLE, coords_data);

  FmsComponentSetCoordinates(volume, coords);
}


/* -------------------------------------------------------------------------- */
/* Construct an FMS data collection with two hexahedra */
/* -------------------------------------------------------------------------- */

void TestMakeFmsHexMesh(FmsDataCollection *dc_ptr) {
  FmsMesh mesh;
  FmsMeshConstruct(&mesh);

  FmsDomain *domains;
  FmsMeshAddDomains(mesh, "domain", 1, &domains);
  FmsDomainSetNumVertices(domains[0], 12);

  FmsDomainSetNumEntities(domains[0], FMS_EDGE, FMS_INT32, 20);
  const int edge_vert[] = {
    0,1,  // 1
    1,2,  // 2
    0,3,  // 3
    1,4,  // 4
    2,5,  // 5
    3,4,  // 6
    4,5,  // 7
    0,6,  // 8
    1,7,  // 9
    2,8,  // 10
    3,9,  // 11
    4,10, // 12
    5,11, // 13
    6,7,  // 14
    7,8,  // 15
    6,9,  // 16
    7,10, // 17
    8,11, // 18
    9,10, // 19
    10,11 // 20
  };
  FmsDomainAddEntities(domains[0], FMS_EDGE, NULL, FMS_INT32, edge_vert, 20);

  FmsDomainSetNumEntities(domains[0], FMS_QUADRILATERAL, FMS_INT32, 11);
  const int face_edge[] = {
    0,3,5,2,     // 1
    1,4,6,3,     // 2
    7,2,10,15,   // 3
    8,3,11,16,   // 4
    9,17,12,4,   // 5
    13,15,18,16, // 6
    14,16,19,17, // 7
    0,8,13,7,    // 8
    1,9,14,8,    // 9
    5,11,18,10,  // 10
    6,12,19,11   // 11
  };
  FmsDomainAddEntities(domains[0], FMS_QUADRILATERAL, NULL, FMS_INT32,
                       face_edge, 11);

  FmsDomainSetNumEntities(domains[0], FMS_HEXAHEDRON, FMS_INT32, 2);
  const int hex_face[] = {
    7,9,0,5,2,3,
    8,10,1,6,3,4,
  };
  FmsDomainAddEntities(domains[0], FMS_HEXAHEDRON, NULL, FMS_INT32, hex_face,
                       2);

  FmsComponent volume, boundary;
  FmsMeshAddComponent(mesh, "volume", &volume);
  FmsComponentAddDomain(volume, domains[0]);

  FmsMeshAddComponent(mesh, "boundary", &boundary);
  FmsInt part_id;
  FmsComponentAddPart(boundary, domains[0], &part_id);
  const int bdr_face[] = {0,1,2,4,5,6,7,8,9,10};
  FmsComponentAddPartEntities(boundary, part_id, FMS_QUADRILATERAL, FMS_INT32,
                              FMS_INT32, FMS_INT32, NULL, bdr_face, NULL, 10);

  FmsComponentAddRelation(volume, 1); // 1 = index of "boundary" component

  FmsTag attributes, bdr_attributes;
  FmsMeshAddTag(mesh, "attributes", &attributes);
  FmsTagSetComponent(attributes, volume);
  const int hex_attr[] = {1,2};
  FmsTagSet(attributes, FMS_UINT8, FMS_INT32, hex_attr, 2);
  const int described_attr[] = {1,2};
  const char *attr_descr[] = {"material 1", "material 2"};
  FmsTagAddDescriptions(attributes, FMS_INT32, described_attr, attr_descr, 2);

  FmsMeshAddTag(mesh, "bdr_attributes", &bdr_attributes);
  FmsTagSetComponent(bdr_attributes, boundary);
  const int bdr_attr[] = {1,1,2,3,4,4,5,5,6,6};
  FmsTagSet(bdr_attributes, FMS_UINT8, FMS_INT32, bdr_attr, 10);
  const int bdr_described_attr[] = {1,2,3,4,5,6};
  const char *bdr_attr_descr[] = {
    "face 1", "face 2", "face 3", "face 4", "face 5", "face 6"
  };
  FmsTagAddDescriptions(bdr_attributes, FMS_INT32, bdr_described_attr,
                        bdr_attr_descr, 6);

  FmsMeshFinalize(mesh);
  FmsMeshValidate(mesh);

  FmsDataCollectionCreate(mesh, "data collection", dc_ptr);
  FmsDataCollection dc = *dc_ptr;

  FmsFieldDescriptor coords_fd;
  FmsDataCollectionAddFieldDescriptor(dc, "coords descriptor", &coords_fd);
  FmsFieldDescriptorSetComponent(coords_fd, volume);
  FmsInt order = 1;
  FmsFieldDescriptorSetFixedOrder(coords_fd, FMS_CONTINUOUS,
                                  FMS_NODAL_GAUSS_CLOSED, order);

  FmsField coords;
  FmsDataCollectionAddField(dc, "coords", &coords);
  const double coords_data[] = {
    0.,0.,0.,
    1.,0.,0.,
    2.,0.,0.,
    0.,0.,1.,
    1.,0.,1.,
    2.,0.,1.,
    0.,1.,0.,
    1.,1.,0.,
    2.,1.,0.,
    0.,1.,1.,
    1.,1.,1.,
    2.,1.,1.
  };
  FmsFieldSet(coords, coords_fd, 3, FMS_BY_VDIM, FMS_DOUBLE, coords_data);

  FmsComponentSetCoordinates(volume, coords);
}


/* -------------------------------------------------------------------------- */
/* FMS to MFEM conversion and testing function */
/* -------------------------------------------------------------------------- */

void TestFmsToMfem(FmsDataCollection dc) {
  mfem::Mesh *mesh;
  int err = ConvertFmsToMfem(dc, &mesh);

  if (err) {
    std::cout << "\nError in ConvertFmsToMfem(), error code: " << err << "\n\n";
    return;
  }

  std::cout << "FMS to MFEM:\n"
            "  Number of vertices     : " << mesh->GetNV() << "\n"
            "  Number of elements     : " << mesh->GetNE() << "\n"
            "  Number of bdr elements : " << mesh->GetNBE() << "\n"
            "  Mesh curavture type    : "
            << mesh->GetNodes()->FESpace()->FEColl()->Name() << "\n"
            "  Number of mesh nodes   : "
            << mesh->GetNodes()->FESpace()->GetNDofs() << '\n';

  // Visualize the mesh
  {
    char vishost[] = "localhost";
    int  visport   = 19916;
    std::cout << "Sending the mesh to GLVis at " << vishost << ':' << visport
              << " ...";
    mfem::socketstream sol_sock(vishost, visport);
    if (sol_sock.good()) {
      sol_sock.precision(8);
      sol_sock << "mesh\n" << *mesh << std::flush;
      std::cout << " done.\n";
    } else {
      std::cout << " unable to connect.\n";
    }
  }

  delete mesh;
}


/* -------------------------------------------------------------------------- */
/* FMS to MFEM conversion function */
/* -------------------------------------------------------------------------- */

int ConvertFmsToMfem(FmsDataCollection dc, mfem::Mesh **mesh_p) {

  FmsInt dim, n_vert, n_elem, n_bdr_elem, space_dim;

  FmsMesh fms_mesh;
  FmsDataCollectionGetMesh(dc, &fms_mesh);

  // Find the first component that has coordinates - that will be the new mfem
  // mesh.
  FmsInt num_comp;
  FmsMeshGetNumComponents(fms_mesh, &num_comp);
  FmsComponent main_comp = NULL;
  FmsField coords = NULL;
  for (FmsInt comp_id = 0; comp_id < num_comp; comp_id++) {
    FmsComponent comp;
    FmsMeshGetComponent(fms_mesh, comp_id, &comp);
    FmsComponentGetCoordinates(comp, &coords);
    if (coords) { main_comp = comp; break; }
  }
  if (!main_comp) { return 1; }
  FmsComponentGetDimension(main_comp, &dim);
  FmsComponentGetNumEntities(main_comp, &n_elem);
  FmsInt n_ents[FMS_NUM_ENTITY_TYPES];
  FmsInt n_main_parts;
  FmsComponentGetNumParts(main_comp, &n_main_parts);
  for (FmsInt et = FMS_VERTEX; et < FMS_NUM_ENTITY_TYPES; et++) {
    n_ents[et] = 0;
    for (FmsInt part_id = 0; part_id < n_main_parts; part_id++) {
      FmsInt num_ents;
      FmsComponentGetPart(main_comp, part_id, (FmsEntityType)et, NULL, NULL,
                          NULL, NULL, &num_ents);
      n_ents[et] += num_ents;
    }
  }
  n_vert = n_ents[FMS_VERTEX];

  // The first related component of dimension dim-1 will be the boundary of the
  // new mfem mesh.
  FmsComponent bdr_comp = NULL;
  FmsInt num_rel_comps;
  const FmsInt *rel_comp_ids;
  FmsComponentGetRelations(main_comp, &rel_comp_ids, &num_rel_comps);
  for (FmsInt i = 0; i < num_rel_comps; i++) {
    FmsComponent comp;
    FmsMeshGetComponent(fms_mesh, rel_comp_ids[i], &comp);
    FmsInt comp_dim;
    FmsComponentGetDimension(comp, &comp_dim);
    if (comp_dim == dim-1) { bdr_comp = comp; break; }
  }
  if (bdr_comp) {
    FmsComponentGetNumEntities(bdr_comp, &n_bdr_elem);
  } else {
    n_bdr_elem = 0;
  }

  FmsFieldGet(coords, NULL, &space_dim, NULL, NULL, NULL);

  int err = 0;
  mfem::Mesh *mesh;
  mesh = new mfem::Mesh(dim, n_vert, n_elem, n_bdr_elem, space_dim);

  // Element tags
  FmsInt num_tags;
  FmsMeshGetNumTags(fms_mesh, &num_tags);
  FmsTag elem_tag = NULL, bdr_tag = NULL;
  for (FmsInt tag_id = 0; tag_id < num_tags; tag_id++) {
    FmsTag tag;
    FmsMeshGetTag(fms_mesh, tag_id, &tag);
    FmsComponent comp;
    FmsTagGetComponent(tag, &comp);
    if (!elem_tag && comp == main_comp) {
      elem_tag = tag;
    } else if (!bdr_tag && comp == bdr_comp) {
      bdr_tag = tag;
    }
  }
  FmsIntType attr_type;
  const void *v_attr, *v_bdr_attr;
  mfem::Array<int> attr, bdr_attr;
  FmsInt num_attr;
  // Element attributes
  if (elem_tag) {
    FmsTagGet(elem_tag, &attr_type, &v_attr, &num_attr);
    if (attr_type == FMS_UINT8) {
      mfem::Array<uint8_t> at((uint8_t*)v_attr, num_attr);
      attr = at;
    } else if (attr_type == FMS_INT32 || attr_type == FMS_UINT32) {
      attr.MakeRef((int*)v_attr, num_attr);
    } else {
      err = 1; // "attribute type not supported!"
      goto func_exit;
    }
  }
  // Boundary attributes
  if (bdr_tag) {
    FmsTagGet(bdr_tag, &attr_type, &v_bdr_attr, &num_attr);
    if (attr_type == FMS_UINT8) {
      mfem::Array<uint8_t> at((uint8_t*)v_bdr_attr, num_attr);
      bdr_attr = at;
    } else if (attr_type == FMS_INT32 || attr_type == FMS_UINT32) {
      bdr_attr.MakeRef((int*)v_bdr_attr, num_attr);
    } else {
      err = 2; // "bdr attribute type not supported!"
      goto func_exit;
    }
  }

  // Add elements
  for (FmsInt part_id = 0; part_id < n_main_parts; part_id++) {
    for (int et = FMS_VERTEX; et < FMS_NUM_ENTITY_TYPES; et++) {
      if (FmsEntityDim[et] != dim) { continue; }

      FmsDomain domain;
      FmsIntType elem_id_type;
      const void *elem_ids;
      const FmsOrientation *elem_ori;
      FmsInt num_elems;
      FmsComponentGetPart(main_comp, part_id, (FmsEntityType)et, &domain,
                          &elem_id_type, &elem_ids, &elem_ori, &num_elems);
      if (num_elems == 0) { continue; }

      if (elem_ids != NULL &&
          (elem_id_type != FMS_INT32 && elem_id_type != FMS_UINT32)) {
        err = 3; goto func_exit;
      }
      if (elem_ori != NULL) {
        err = 4; goto func_exit;
      }

      const FmsInt nv = FmsEntityNumVerts[et];
      mfem::Array<int> ents_verts(num_elems*nv);
      if (elem_ids == NULL) {
        FmsDomainGetEntitiesVerts(domain, (FmsEntityType)et, NULL, FMS_INT32,
                                  0, ents_verts.GetData(), num_elems);
      } else {
        const int *ei = (const int *)elem_ids;
        for (FmsInt i = 0; i < num_elems; i++) {
          FmsDomainGetEntitiesVerts(domain, (FmsEntityType)et, NULL, FMS_INT32,
                                    ei[i], &ents_verts[i*nv], 1);
        }
      }
      const int elem_offset = mesh->GetNE();
      switch ((FmsEntityType)et) {
      case FMS_EDGE:
        err = 5;
        goto func_exit;
        break;
      case FMS_TRIANGLE:
        for (FmsInt i = 0; i < num_elems; i++) {
          mesh->AddTriangle(
            &ents_verts[3*i], elem_tag ? attr[elem_offset+i] : 1);
        }
        break;
      case FMS_QUADRILATERAL:
        for (FmsInt i = 0; i < num_elems; i++) {
          mesh->AddQuad(&ents_verts[4*i], elem_tag ? attr[elem_offset+i] : 1);
        }
        break;
      case FMS_TETRAHEDRON:
        for (FmsInt i = 0; i < num_elems; i++) {
          mesh->AddTet(&ents_verts[4*i], elem_tag ? attr[elem_offset+i] : 1);
        }
        break;
      case FMS_HEXAHEDRON:
        for (FmsInt i = 0; i < num_elems; i++) {
          mesh->AddHex(&ents_verts[8*i], elem_tag ? attr[elem_offset+i] : 1);
        }
        break;
      default:
        break;
      }
    }
  }

  // Add boundary elements
  if (bdr_comp && n_bdr_elem > 0) {
    FmsInt n_bdr_parts;
    FmsComponentGetNumParts(bdr_comp, &n_bdr_parts);

    for (FmsInt part_id = 0; part_id < n_bdr_parts; part_id++) {
      for (int et = FMS_VERTEX; et < FMS_NUM_ENTITY_TYPES; et++) {
        if (FmsEntityDim[et] != dim-1) { continue; }

        FmsDomain domain;
        FmsIntType elem_id_type;
        const void *elem_ids;
        const FmsOrientation *elem_ori;
        FmsInt num_elems;
        FmsComponentGetPart(bdr_comp, part_id, (FmsEntityType)et, &domain,
                            &elem_id_type, &elem_ids, &elem_ori, &num_elems);
        if (num_elems == 0) { continue; }

        if (elem_ids != NULL &&
            (elem_id_type != FMS_INT32 && elem_id_type != FMS_UINT32)) {
          err = 6; goto func_exit;
        }
        if (elem_ori != NULL) {
          err = 7; goto func_exit;
        }

        const FmsInt nv = FmsEntityNumVerts[et];
        mfem::Array<int> ents_verts(num_elems*nv);
        if (elem_ids == NULL) {
          FmsDomainGetEntitiesVerts(domain, (FmsEntityType)et, NULL, FMS_INT32,
                                    0, ents_verts.GetData(), num_elems);
        } else {
          const int *ei = (const int *)elem_ids;
          for (FmsInt i = 0; i < num_elems; i++) {
            FmsDomainGetEntitiesVerts(domain, (FmsEntityType)et, NULL,
                                      FMS_INT32, ei[i], &ents_verts[i*nv], 1);
          }
        }
        const int elem_offset = mesh->GetNBE();
        switch ((FmsEntityType)et) {
        case FMS_EDGE:
          for (FmsInt i = 0; i < num_elems; i++) {
            mesh->AddBdrSegment(
              &ents_verts[2*i], bdr_tag ? bdr_attr[elem_offset+i] : 1);
          }
          break;
        case FMS_TRIANGLE:
          for (FmsInt i = 0; i < num_elems; i++) {
            mesh->AddBdrTriangle(
              &ents_verts[3*i], bdr_tag ? bdr_attr[elem_offset+i] : 1);
          }
          break;
        case FMS_QUADRILATERAL:
          for (FmsInt i = 0; i < num_elems; i++) {
            mesh->AddBdrQuad(
              &ents_verts[4*i], bdr_tag ? bdr_attr[elem_offset+i] : 1);
          }
          break;
        default:
          break;
        }
      }
    }
  }

  // Finalize the mesh topology
  // FIXME: mfem::Mesh::FinalizeCheck() assumes all vertices are added
  // mesh.FinalizeTopology();

  // Transfer coordinates
  {
    FmsFieldDescriptor coords_fd;
    FmsLayoutType coords_layout;
    FmsScalarType coords_data_type;
    const void *coords_data;
    FmsFieldGet(coords, &coords_fd, NULL, &coords_layout, &coords_data_type,
                &coords_data);

    FmsInt coords_num_dofs;
    FmsFieldDescriptorGetNumDofs(coords_fd, &coords_num_dofs);

    // assuming coordinates data type is double
    if (coords_data_type != FMS_DOUBLE) {
      err = 8; goto func_exit;
    }
    const double *coords_dbl_data = (const double *)coords_data;

    FmsFieldDescriptorType coords_fd_type;
    FmsFieldDescriptorGetType(coords_fd, &coords_fd_type);
    if (coords_fd_type != FMS_FIXED_ORDER) {
      err = 9; goto func_exit;
    }
    FmsFieldType coords_field_type;
    FmsBasisType coords_basis_type;
    FmsInt coords_order;
    FmsFieldDescriptorGetFixedOrder(coords_fd, &coords_field_type,
                                    &coords_basis_type, &coords_order);
    if (coords_field_type != FMS_CONTINUOUS) {
      err = 10; goto func_exit;
    }
    if (coords_basis_type != FMS_NODAL_GAUSS_CLOSED) {
      err = 11; goto func_exit;
    }

    const FmsInt nstride = (coords_layout == FMS_BY_VDIM) ? space_dim : 1;
    const FmsInt vstride = (coords_layout == FMS_BY_VDIM) ? 1 : coords_num_dofs;

    // Set the vertex coordinates to zero
    const double origin[3] = {0.,0.,0.};
    for (FmsInt vi = 0; vi < n_vert; vi++) {
      mesh->AddVertex(origin);
    }

    // Finalize the mesh topology
    mesh->FinalizeTopology();

    // Switch to mfem::Mesh with nodes (interpolates the linear coordinates)
    const bool discont = false;
    mesh->SetCurvature(coords_order, discont, space_dim,
                       (coords_layout == FMS_BY_VDIM) ?
                       mfem::Ordering::byVDIM : mfem::Ordering::byNODES);

    // Finalize mesh construction
    mesh->Finalize();

    // Set the high-order mesh nodes
    mfem::GridFunction &nodes = *mesh->GetNodes();
    if ((FmsInt)(nodes.Size()) != coords_num_dofs*space_dim) {
      err = 12; goto func_exit;
    }
    mfem::FiniteElementSpace *fes = nodes.FESpace();
    const int vdim = fes->GetVDim();
    const mfem::FiniteElementCollection *fec = fes->FEColl();
    const int vert_dofs = fec->DofForGeometry(mfem::Geometry::POINT);
    const int edge_dofs = fec->DofForGeometry(mfem::Geometry::SEGMENT);
    const int tri_dofs = fec->DofForGeometry(mfem::Geometry::TRIANGLE);
    const int quad_dofs = fec->DofForGeometry(mfem::Geometry::SQUARE);
    const int tet_dofs = fec->DofForGeometry(mfem::Geometry::TETRAHEDRON);
    const int hex_dofs = fec->DofForGeometry(mfem::Geometry::CUBE);
    int ent_dofs[FMS_NUM_ENTITY_TYPES];
    ent_dofs[FMS_VERTEX] = vert_dofs;
    ent_dofs[FMS_EDGE] = edge_dofs;
    ent_dofs[FMS_TRIANGLE] = tri_dofs;
    ent_dofs[FMS_QUADRILATERAL] = quad_dofs;
    ent_dofs[FMS_TETRAHEDRON] = tet_dofs;
    ent_dofs[FMS_HEXAHEDRON] = hex_dofs;
    FmsInt fms_dof_offset = 0;
    int mfem_ent_cnt[4] = {0,0,0,0}; // mfem entity counters, by dimension
    int mfem_last_vert_cnt = 0;
    mfem::HashTable<mfem::Hashed2> mfem_edge;
    mfem::HashTable<mfem::Hashed4> mfem_face;
    if (dim >= 2 && edge_dofs > 0) {
      mfem::Array<int> ev;
      for (int i = 0; i < mesh->GetNEdges(); i++) {
        mesh->GetEdgeVertices(i, ev);
        int id = mfem_edge.GetId(ev[0], ev[1]);
        if (id != i) { err = 13; goto func_exit; }
      }
    }
    if (dim >= 3 &&
        ((n_ents[FMS_TRIANGLE] > 0 && tri_dofs > 0) ||
         (n_ents[FMS_QUADRILATERAL] > 0 && quad_dofs > 0))) {
      mfem::Array<int> fv;
      for (int i = 0; i < mesh->GetNFaces(); i++) {
        mesh->GetFaceVertices(i, fv);
        if (fv.Size() == 3) { fv.Append(INT_MAX); }
        // HashTable uses the smallest 3 of the 4 indices to hash Hashed4
        int id = mfem_face.GetId(fv[0], fv[1], fv[2], fv[3]);
        if (id != i) { err = 14; goto func_exit; }
      }
    }

    // Loop over all parts of the main component
    for (FmsInt part_id = 0; part_id < n_main_parts; part_id++) {
      // Loop over all entity types in the part
      for (FmsInt et = FMS_VERTEX; et < FMS_NUM_ENTITY_TYPES; et++) {
        FmsDomain domain;
        FmsIntType ent_id_type;
        const void *ents;
        const FmsOrientation *ents_ori;
        FmsInt num_ents;
        FmsComponentGetPart(main_comp, part_id, (FmsEntityType)et, &domain,
                            &ent_id_type, &ents, &ents_ori, &num_ents);
        if (num_ents == 0) { continue; }
        if (ent_dofs[et] == 0) {
          if (et == FMS_VERTEX) { mfem_last_vert_cnt = mfem_ent_cnt[et]; }
          mfem_ent_cnt[FmsEntityDim[et]] += num_ents;
          continue;
        }

        if (ents != NULL &&
            (ent_id_type != FMS_INT32 && ent_id_type != FMS_UINT32)) {
          err = 15; goto func_exit;
        }
        if (ents_ori != NULL) {
          err = 16; goto func_exit;
        }

        if (et == FMS_VERTEX) {
          const int mfem_dof_offset = mfem_ent_cnt[0]*vert_dofs;
          for (FmsInt i = 0; i < num_ents*vert_dofs; i++) {
            for (int j = 0; j < vdim; j++) {
              const int idx = i*nstride+j*vstride;
              nodes(mfem_dof_offset*nstride+idx) =
                coords_dbl_data[fms_dof_offset*nstride+idx];
            }
          }
          fms_dof_offset += num_ents*vert_dofs;
          mfem_last_vert_cnt = mfem_ent_cnt[et];
          mfem_ent_cnt[0] += num_ents;
          continue;
        }
        mfem::Array<int> dofs;
        if (FmsEntityDim[et] == dim) {
          for (FmsInt e = 0; e < num_ents; e++) {
            fes->GetElementInteriorDofs(mfem_ent_cnt[dim]+e, dofs);
            for (int i = 0; i < ent_dofs[et]; i++, fms_dof_offset++) {
              for (int j = 0; j < vdim; j++) {
                nodes(fes->DofToVDof(dofs[i],j)) =
                  coords_dbl_data[fms_dof_offset*nstride+j*vstride];
              }
            }
          }
          mfem_ent_cnt[dim] += num_ents;
          continue;
        }
        const FmsInt nv = FmsEntityNumVerts[et];
        mfem::Array<int> ents_verts(num_ents*nv), m_ev;
        const int *ei = (const int *)ents;
        if (ents == NULL) {
          FmsDomainGetEntitiesVerts(domain, (FmsEntityType)et, NULL, FMS_INT32,
                                    0, ents_verts.GetData(), num_ents);
        } else {
          for (FmsInt i = 0; i < num_ents; i++) {
            FmsDomainGetEntitiesVerts(domain, (FmsEntityType)et, NULL,
                                      FMS_INT32, ei[i], &ents_verts[i*nv], 1);
          }
        }
        for (int i = 0; i < ents_verts.Size(); i++) {
          ents_verts[i] += mfem_last_vert_cnt;
        }
        const int *perm;
        switch ((FmsEntityType)et) {
        case FMS_EDGE: {
          for (FmsInt part_ent_id = 0; part_ent_id < num_ents; part_ent_id++) {
            const int *ev = &ents_verts[2*part_ent_id];
            int mfem_edge_id = mfem_edge.FindId(ev[0], ev[1]);
            if (mfem_edge_id < 0) {
              err = 17; goto func_exit;
            }
            mesh->GetEdgeVertices(mfem_edge_id, m_ev);
            int ori = (ev[0] == m_ev[0]) ? 0 : 1;
            perm = fec->DofOrderForOrientation(mfem::Geometry::SEGMENT, ori);
            fes->GetEdgeInteriorDofs(mfem_edge_id, dofs);
            for (int i = 0; i < edge_dofs; i++) {
              for (int j = 0; j < vdim; j++) {
                nodes(fes->DofToVDof(dofs[i],j)) =
                  coords_dbl_data[(fms_dof_offset+perm[i])*nstride+j*vstride];
              }
            }
            fms_dof_offset += edge_dofs;
          }
          break;
        }
        case FMS_TRIANGLE: {
          for (FmsInt part_ent_id = 0; part_ent_id < num_ents; part_ent_id++) {
            const int *tv = &ents_verts[3*part_ent_id];
            int mfem_face_id = mfem_face.FindId(tv[0], tv[1], tv[2], INT_MAX);
            if (mfem_face_id < 0) {
              err = 18; goto func_exit;
            }
            mesh->GetFaceVertices(mfem_face_id, m_ev);
            int ori = 0;
            while (tv[ori] != m_ev[0]) { ori++; }
            ori = (tv[(ori+1)%3] == m_ev[1]) ? 2*ori : 2*ori+1;
            perm = fec->DofOrderForOrientation(mfem::Geometry::TRIANGLE, ori);
            fes->GetFaceInteriorDofs(mfem_face_id, dofs);
            for (int i = 0; i < tri_dofs; i++) {
              for (int j = 0; j < vdim; j++) {
                nodes(fes->DofToVDof(dofs[i],j)) =
                  coords_dbl_data[(fms_dof_offset+perm[i])*nstride+j*vstride];
              }
            }
            fms_dof_offset += tri_dofs;
          }
          break;
        }
        case FMS_QUADRILATERAL: {
          for (FmsInt part_ent_id = 0; part_ent_id < num_ents; part_ent_id++) {
            const int *qv = &ents_verts[4*part_ent_id];
            int mfem_face_id = mfem_face.FindId(qv[0], qv[1], qv[2], qv[3]);
            if (mfem_face_id < 0) { err = 19; goto func_exit; }
            mesh->GetFaceVertices(mfem_face_id, m_ev);
            int ori = 0;
            while (qv[ori] != m_ev[0]) { ori++; }
            ori = (qv[(ori+1)%4] == m_ev[1]) ? 2*ori : 2*ori+1;
            perm = fec->DofOrderForOrientation(mfem::Geometry::SQUARE, ori);
            fes->GetFaceInteriorDofs(mfem_face_id, dofs);
            for (int i = 0; i < quad_dofs; i++) {
              for (int j = 0; j < vdim; j++) {
                nodes(fes->DofToVDof(dofs[i],j)) =
                  coords_dbl_data[(fms_dof_offset+perm[i])*nstride+j*vstride];
              }
            }
            fms_dof_offset += quad_dofs;
          }
          break;
        }
        default: break;
        }
        mfem_ent_cnt[FmsEntityDim[et]] += num_ents;
      }
    }
  }

func_exit:

  if (err) {
    delete mesh;
  } else {
    *mesh_p = mesh;
  }
  return err;
}
