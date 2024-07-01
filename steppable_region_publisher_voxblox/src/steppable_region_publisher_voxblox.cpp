#include <steppable_region_publisher_voxblox/steppable_region_publisher_voxblox.h>


namespace steppable_region_publisher_voxblox {

  bool calcSteppableRegion(const std::shared_ptr<voxblox::EsdfMap>& esdf_map,
                           const std::shared_ptr<voxblox::MeshLayer>& mesh_layer,
                           std::shared_ptr<voxblox::Mesh>& steppable_region,
                           const SurfaceParam& param
                           ) {
    if(!steppable_region){
      steppable_region = std::make_shared<voxblox::Mesh>();
    }
    steppable_region->clear();

    if(!esdf_map || !mesh_layer){
      std::cerr << __FUNCTION__ << " invalid input!" << std::endl;
      return false;
    }

    std::vector<Eigen::Vector3d> supportSole;
    std::vector<Eigen::Vector3d> collisionSole;
    for(double x=0.0; x<=std::max(param.collisionSoleRange,param.supportSoleRange); x+=param.soleResolution){
      for(double y=0.0; y<=std::max(param.collisionSoleRange,param.supportSoleRange); y+=param.soleResolution){
        if(std::pow(x,2) + std::pow(y,2) <= std::pow(param.supportSoleRange,2)){
          supportSole.emplace_back(x,y,0.0);
        }else if(std::pow(x,2) + std::pow(y,2) <= std::pow(param.collisionSoleRange,2)){
          collisionSole.emplace_back(x,y,0.0);
        }
      }
    }

    voxblox::BlockIndexList mesh_indices;
    mesh_layer->getAllAllocatedMeshes(&mesh_indices);
    int idx = 0;
    for(int m=0; m<mesh_indices.size(); m++){
      std::shared_ptr<const voxblox::Mesh> mesh = mesh_layer->getMeshPtrByIndex(mesh_indices[m]);

      for(int f=0;f+2<mesh->indices.size();f+=3){
        Eigen::Vector3f n0 = mesh->normals[mesh->indices[f]];
        Eigen::Vector3f n1 = mesh->normals[mesh->indices[f+1]];
        Eigen::Vector3f n2 = mesh->normals[mesh->indices[f+2]];
        Eigen::Vector3d n = ((n0 + n1 + n2) / 3).normalized().cast<double>();
        if(param.maxNormalAngle < std::acos(std::min(1.0,(std::max(-1.0,n.dot(Eigen::Vector3d::UnitZ())))))){
          continue;
        }

        Eigen::Vector3f p0 = mesh->vertices[mesh->indices[f]];
        Eigen::Vector3f p1 = mesh->vertices[mesh->indices[f+1]];
        Eigen::Vector3f p2 = mesh->vertices[mesh->indices[f+2]];
        Eigen::Vector3d p = ((p0 + p1 + p2) / 3).cast<double>();

        Eigen::Matrix3d R;
        {
          Eigen::Vector3d x,y,z;
          z = n;
          if(Eigen::Vector3d::UnitY().cross(z).norm() > 0){
            x = Eigen::Vector3d::UnitY().cross(z).normalized();
            y = z.cross(x).normalized();
          }else{
            y = z.cross(Eigen::Vector3d::UnitX()).normalized();
            x = y.cross(z).normalized();
          }
          R.col(0) = x; R.col(1) = y; R.col(2) = z;
        }

        bool steppable = false;
        for(double offset = 0.0; !steppable && offset <= param.maxSupportDistance; offset+=param.offsetResolution){
          Eigen::Vector3d origin = p + offset * n;
          bool valid = true;
          for(int s = 0; valid && s<supportSole.size(); s++){
            Eigen::Vector3d sample = origin + R * supportSole[s];
            double d;
            if(!esdf_map->getDistanceAtPosition(sample, &d) ||
               d > param.maxSupportDistance ||
               d < param.minCollisionDistance) {
              valid = false;
            }
          }
          for(int s = 0; valid && s<collisionSole.size(); s++){
            Eigen::Vector3d sample = origin + R * collisionSole[s];
            double d;
            if(!esdf_map->getDistanceAtPosition(sample, &d) ||
               d < param.minCollisionDistance) {
              valid = false;
            }
          }
          if(valid) steppable = true;
        }
        if(!steppable) continue;

        steppable_region->vertices.push_back(p0);
        steppable_region->vertices.push_back(p1);
        steppable_region->vertices.push_back(p2);
        steppable_region->indices.push_back(idx);
        steppable_region->indices.push_back(idx+1);
        steppable_region->indices.push_back(idx+2);
        steppable_region->normals.push_back(n0);
        steppable_region->normals.push_back(n1);
        steppable_region->normals.push_back(n2);
        idx += 3;
      }
    }
    return true;
  }

};
