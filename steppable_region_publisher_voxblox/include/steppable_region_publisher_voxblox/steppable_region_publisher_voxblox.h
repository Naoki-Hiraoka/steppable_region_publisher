#ifndef STEPPABLE_REGION_PUBLISHER_VOXBLOX_STEPPABLE_REGION_PUBLISHER_VOXBLOX_H
#define STEPPABLE_REGION_PUBLISHER_VOXBLOX_STEPPABLE_REGION_PUBLISHER_VOXBLOX_H

#include <voxblox/core/esdf_map.h>
#include <voxblox/mesh/mesh_layer.h>

namespace steppable_region_publisher_voxblox {

  class SurfaceParam {
  public:
    double collisionSoleRange = 0.15; // endeffector相対半径
    double minCollisionDistance = 0.0;
    double supportSoleRange = 0.1; // endeffector相対半径
    double maxSupportDistance = 0.02;
    double soleResolution = 0.02;

    float offsetResolution = 0.005;

    double maxNormalAngle = 30.0 / 180.0 * M_PI;
  };

  bool calcSteppableRegion(const std::shared_ptr<voxblox::EsdfMap>& esdf_map,
                           const std::shared_ptr<voxblox::MeshLayer>& mesh_layer,
                           std::shared_ptr<voxblox::Mesh>& steppable_region,
                           const SurfaceParam& param = SurfaceParam()
                           );

};

#endif
