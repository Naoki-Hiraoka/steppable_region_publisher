#include <choreonoid_viewer/choreonoid_viewer.h>
#include <cnoid/Body>
#include <cnoid/BodyLoader>
#include <cnoid/SceneMarkers>
#include <iostream>
#include <ros/package.h>

#include <cnoid/TimeMeasure>
#include <cnoid/MeshExtractor>
#include <cnoid/MeshFilter>

#include <steppable_region_publisher_voxblox/steppable_region_publisher_voxblox.h>
#include <choreonoid_voxblox/choreonoid_voxblox.h>

#include <voxblox/core/tsdf_map.h>
#include <voxblox/core/esdf_map.h>
#include <voxblox/integrator/tsdf_integrator.h>
#include <voxblox/integrator/esdf_integrator.h>
#include <voxblox/mesh/mesh_integrator.h>
#include <voxblox/utils/planning_utils.h>
#include <voxblox/io/layer_io.h>

namespace steppable_region_publisher_voxblox{
  inline std::vector<std::pair<cnoid::Vector3, cnoid::Vector3> > getSurfaceVerticesAndNormals(cnoid::LinkPtr link, float resolution, float minangle) {
    // 1つのvertexを取得したら、resolutionのサイズの同じ立方体の中にありかつ法線がminangle以下の他のvertexは取得しない
    // faceが巨大な場合、faceの内部の点をresolutionの間隔でサンプリングして取得する

    cnoid::MeshExtractor meshExtractor;
    cnoid::MeshFilter meshFilter;

    std::vector<std::pair<cnoid::Vector3, cnoid::Vector3> > vertices;
    cnoid::SgMeshPtr mesh = meshExtractor.integrate(link->collisionShape());
    if(mesh && (mesh->numTriangles() != 0)) {
      meshFilter.generateNormals(mesh,M_PI,true);
      mesh->updateBoundingBox();
      cnoid::BoundingBoxf bbx = mesh->boundingBox();
      cnoid::Vector3f bbxSize = bbx.max() - bbx.min();
      std::vector<std::vector<std::vector<std::vector<cnoid::Vector3f> > > > bin; // normalを入れる
      bin.resize(int(bbxSize[0]/resolution)+1);
      for(int x=0;x<bin.size();x++){
        bin[x].resize(int(bbxSize[1]/resolution)+1);
        for(int y=0;y<bin[x].size();y++){
          bin[x][y].resize(int(bbxSize[2]/resolution)+1);
        }
      }

      for(int j=0;j<mesh->numTriangles();j++){
        cnoid::Vector3f v0 = mesh->vertices()->at(mesh->triangle(j)[0]);
        cnoid::Vector3f v1 = mesh->vertices()->at(mesh->triangle(j)[1]);
        cnoid::Vector3f v2 = mesh->vertices()->at(mesh->triangle(j)[2]);
        cnoid::Vector3f normal = mesh->normals()->at(mesh->normalIndices()[j*3]); // linkの外側に向かう方向
        float l1 = (v1 - v0).norm();
        float l2 = (v2 - v0).norm();
        cnoid::Vector3f d1 = (v1 - v0).normalized();
        cnoid::Vector3f d2 = (v2 - v0).normalized();

        for(float m=0;;){
          float n_max = (l1==0)? l2 : l2*(1-m/l1);
          for(float n=0;;){
            cnoid::Vector3f v = v0 + d1 * m + d2 * n;
            int x = int((v[0] - bbx.min()[0])/resolution);
            int y = int((v[1] - bbx.min()[1])/resolution);
            int z = int((v[2] - bbx.min()[2])/resolution);

            bool exists = false;
            for(int s=0;s<bin[x][y][z].size();s++){
              if(minangle >= std::acos(std::min(1.0f,(std::max(-1.0f,normal.dot(bin[x][y][z][s])))))){
                exists = true;
                break;
              }
            }
            if(!exists){
              bin[x][y][z].push_back(normal);
              vertices.emplace_back(v.cast<double>(), normal.cast<double>());
            }

            if(n>= n_max) break;
            else n = std::min(n+resolution, n_max);
          }

          if(m>=l1) break;
          else m = std::min(m+resolution, l1);
        }
      }
    }
    return vertices;
  }

  void sample1_vnoid(){
    cnoid::TimeMeasure timer; timer.begin();

    // load robot
    std::string modelfile = ros::package::getPath("vnoid_world") + "/vnoid/model/field2022/field2022_athletics.body";
    cnoid::BodyLoader bodyLoader;
    cnoid::BodyPtr world = bodyLoader.load(modelfile);
    world->rootLink()->p() = cnoid::Vector3(0,0,0);
    world->calcCenterOfMass();


    // setup viewer
    choreonoid_viewer::Viewer viewer;
    viewer.objects(world);
    viewer.drawObjects();

    double voxel_size = 0.02;
    double radius = 3.0;
    double default_distance = 0.5;
    std::string tsdfFileName = ros::package::getPath("steppable_region_publisher_voxblox") + "/layer/sample1_vnoid_tsdf";
    std::string esdfFileName = ros::package::getPath("steppable_region_publisher_voxblox") + "/layer/sample1_vnoid_esdf";
    voxblox::TsdfMap::Config tsdf_config;
    tsdf_config.tsdf_voxel_size = voxel_size;
    std::shared_ptr<voxblox::Layer<voxblox::TsdfVoxel> > tsdf_layer = std::make_shared<voxblox::Layer<voxblox::TsdfVoxel> >(tsdf_config.tsdf_voxel_size, tsdf_config.tsdf_voxels_per_side);

    if(!voxblox::io::LoadLayer<voxblox::TsdfVoxel>(tsdfFileName,
                                                   &tsdf_layer)){

      // {
      //   voxblox::HierarchicalIndexMap block_voxel_list;
      //   voxblox::utils::getAndAllocateSphereAroundPoint(voxblox::Point(0.0, 0.0, 0.0), radius,
      //                                                   tsdf_layer.get(), &block_voxel_list);
      //   for (const std::pair<voxblox::BlockIndex, voxblox::VoxelIndexList>& kv : block_voxel_list) {
      //     // Get block.
      //     voxblox::Block<voxblox::TsdfVoxel>::Ptr block_ptr = tsdf_layer->getBlockPtrByIndex(kv.first);

      //     for (const voxblox::VoxelIndex& voxel_index : kv.second) {
      //       if (!block_ptr->isValidVoxelIndex(voxel_index)) {
      //         continue;
      //       }
      //       voxblox::TsdfVoxel& tsdf_voxel = block_ptr->getVoxelByVoxelIndex(voxel_index);
      //       tsdf_voxel.distance = default_distance;
      //     }
      //   }
      // }
      voxblox::TsdfIntegratorBase::Config tsdf_integrator_config;
      tsdf_integrator_config.voxel_carving_enabled = true;
      tsdf_integrator_config.default_truncation_distance = voxel_size;
      tsdf_integrator_config.max_weight = 30.0;
      tsdf_integrator_config.min_ray_length_m = 0.01;
      // tsdf_integrator_config.use_sparsity_compensation_factor = true; // simulation用限定. 消えないように
      // tsdf_integrator_config.sparsity_compensation_factor = 5.0;  // simulation用限定
      std::shared_ptr<voxblox::FastTsdfIntegrator> tsdfIntegrator = std::make_shared<voxblox::FastTsdfIntegrator>(tsdf_integrator_config, tsdf_layer.get());
      std::cerr << "after tsdf init: " << timer.measure() << "[s]." << std::endl;

      double ray = 0.1;
      for(int i=0;i<world->numLinks();i++){
        std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d> > vertices = getSurfaceVerticesAndNormals(world->link(i), voxel_size / 2 /*voxsel_sizeよりも多い濃度でsampleしないと穴があきgradientが不正確*/, M_PI/3); // link local
        for(int j=0;j<vertices.size();j++){
          cnoid::Vector3 p = world->link(i)->T() * vertices[j].first; // world frame
          cnoid::Vector3 n = (world->link(i)->R() * vertices[j].second).normalized(); // world frame

          if(p.norm() > radius) continue;

          cnoid::Vector3 origin = p + n * ray; // カメラ原点. world frame
          cnoid::Vector3 z = -n.normalized(); // カメラ姿勢. world frame
          cnoid::Vector3 x, y;
          if(cnoid::Vector3::UnitY().cross(z).norm() > 0){
            x = cnoid::Vector3::UnitY().cross(z).normalized();
            y = z.cross(x).normalized();
          }else{
            y = z.cross(cnoid::Vector3::UnitX()).normalized();
            x = y.cross(z).normalized();
          }
          cnoid::Matrix3 R; R.col(0) = x; R.col(1) = y; R.col(2) = z; // カメラ姿勢. world frame
          voxblox::Transformation trans(Eigen::Quaterniond(R).cast<float>(), origin.cast<float>());
          voxblox::Pointcloud pcl{Eigen::Vector3f(0.0,0.0,ray)};
          voxblox::Colors color{voxblox::Color(0.0,0.0,0.0)};
          tsdfIntegrator->integratePointCloud(trans, pcl, color);
        }
      }

      voxblox::io::SaveLayer(*tsdf_layer,
                             tsdfFileName);
    }
    std::shared_ptr<voxblox::TsdfMap> tsdf_map = std::make_shared<voxblox::TsdfMap>(tsdf_layer);
    std::cerr << "after tsdf: " << timer.measure() << "[s]." << std::endl;

    voxblox::EsdfMap::Config esdf_config;
    esdf_config.esdf_voxel_size = tsdf_config.tsdf_voxel_size;
    std::shared_ptr<voxblox::Layer<voxblox::EsdfVoxel> > esdf_layer = std::make_shared<voxblox::Layer<voxblox::EsdfVoxel> >(esdf_config.esdf_voxel_size, esdf_config.esdf_voxels_per_side);
    if(!voxblox::io::LoadLayer<voxblox::EsdfVoxel>(esdfFileName,
                                                   &esdf_layer)){
      voxblox::EsdfIntegrator::Config esdf_integrator_config;
      esdf_integrator_config.min_distance_m = voxel_size;
      esdf_integrator_config.max_distance_m = default_distance;
      esdf_integrator_config.default_distance_m = default_distance;
      esdf_integrator_config.clear_sphere_radius = radius;
      //esdf_integrator_config.full_euclidean_distance = true;
      std::shared_ptr<voxblox::EsdfIntegrator> esdfIntegrator = std::make_shared<voxblox::EsdfIntegrator>(esdf_integrator_config, tsdf_map->getTsdfLayerPtr(), esdf_layer.get());
      esdfIntegrator->addNewRobotPosition(voxblox::Point(0.0, 0.0, 0.0)); // clearする.
      std::cerr << "after esdf init: " << timer.measure() << "[s]." << std::endl;

      esdfIntegrator->updateFromTsdfLayer(true);

      voxblox::io::SaveLayer(*esdf_layer,
                             esdfFileName);
    }

    std::shared_ptr<voxblox::EsdfMap> esdf_map = std::make_shared<voxblox::EsdfMap>(esdf_layer);
    std::cerr << "after esdf: " << timer.measure() << "[s]." << std::endl;

    std::shared_ptr<voxblox::MeshLayer> mesh_layer = std::make_shared<voxblox::MeshLayer>(tsdf_map->block_size());
    voxblox::MeshIntegratorConfig mesh_config;
    std::shared_ptr<voxblox::MeshIntegrator<voxblox::TsdfVoxel> > meshIntegrator = std::make_shared<voxblox::MeshIntegrator<voxblox::TsdfVoxel> >(mesh_config, tsdf_map->getTsdfLayerPtr(), mesh_layer.get());
    meshIntegrator->generateMesh(true, true);
    std::cerr << "after mesh: " << timer.measure() << "[s]." << std::endl;

    // main loop
    steppable_region_publisher_voxblox::SurfaceParam param;
    param.maxSupportDistance = 0.03;
    param.minCollisionDistance = 0.00;
    param.offsetResolution = 0.01;
    param.collisionSoleRange = 0.161;
    param.supportSoleRange = 0.121;
    param.soleResolution = 0.04;
    std::shared_ptr<voxblox::Mesh> steppable_region;
    steppable_region_publisher_voxblox::calcSteppableRegion(esdf_map, mesh_layer, steppable_region, param);
    std::cerr << "after region: " << timer.measure() << "[s]." << std::endl;

    cnoid::MeshFilter meshFilter;

    cnoid::SgShapePtr outputMesh = new cnoid::SgShape;
    {
      cnoid::SgMeshPtr model = choreonoid_voxblox::convertToChoreonoidModel(steppable_region);
      outputMesh->setMesh(model);
      cnoid::SgMaterialPtr material = new cnoid::SgMaterial();
      material->setTransparency(0);
      cnoid::Vector3f color(255/255.0,128/255.0,255/255.0);
      material->setDiffuseColor(color);
      //material->setSpecularColor(color);
      outputMesh->setMaterial(material);
    }

    cnoid::SgShapePtr worldMesh = new cnoid::SgShape;
    {
      std::shared_ptr<voxblox::Mesh> mesh = std::make_shared<voxblox::Mesh>();
      mesh_layer->getMesh(mesh.get());
      cnoid::SgMeshPtr model = choreonoid_voxblox::convertToChoreonoidModel(mesh);
      worldMesh->setMesh(model);
      cnoid::SgMaterialPtr material = new cnoid::SgMaterial();
      material->setTransparency(0);
      cnoid::Vector3f color(155/255.0,155/255.0,155/255.0);
      material->setDiffuseColor(color);
      //material->setSpecularColor(color);
      worldMesh->setMaterial(material);
    }

    for(int i=0;true;i++){
      std::vector<cnoid::SgNodePtr> markers;
      if(i%2==0){
        markers.push_back(worldMesh);
      }else{
        markers.push_back(outputMesh);
      }
      viewer.drawOn(markers);
      viewer.drawObjects();

      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    return;
  }
}
