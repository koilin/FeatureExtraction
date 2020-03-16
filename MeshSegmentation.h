/*********************************************************************************************************************
* this project is created by Chuhua Xian                                                                             *
* Email: chuhuaxian@gmail.com                                                                                        *
* 2010. 05.                                                                                                          *
**********************************************************************************************************************/

/*********************************************************************************************************************
* File Information:                                                                                                  *
* written by : Chuhua Xian                                                                                           *
* modified by : Chuhua Xian                                                                                          *
* modified date : 2020. 02. 12                                                                                       *
**********************************************************************************************************************/


#ifndef _MESH_SEGMENTATION_H_
#define _MESH_SEGMENTATION_H_

#include <map>
#include <set>
#include <vector>
#include <Core/BaseMeshEntity.h>
#include <Core/PGMeshEntity.h>


namespace GeometryProcess
{
//---------------------------------------------------------------------------------------------------------------------

	class MeshSegmentation
	{
	public:
		// static const std::string STR_MESH_REGION;
		static const std::string STR_SHARP_EDGE;
		typedef std::map<int, std::set<int>> AdjacentGraph;
		typedef std::vector<std::vector<PGMesh::FaceHandle>> MeshRegions;
	public:
		typedef struct SegmentationParameters
		{
			double normal_variance_;
			double threshold_;
			double curvature_variance_;
		}SegmentationParameters;
	public:
		MeshSegmentation();
		~MeshSegmentation();
	public:
		void set_mesh_entity(BaseMeshEntity * _mesh_entity);
		void set_segmentation_parameter(const SegmentationParameters & _param);
		SegmentationParameters get_segmentation_parameter() const;
		int n_regions() const;
		const AdjacentGraph & get_adjacent_graph() const;
		void segment_mesh();
		MeshRegions segment_mesh_entity();

	private:
		BaseMeshEntity * mesh_entity_;
		SegmentationParameters param_;
		AdjacentGraph adjacent_graph_;
		int n_regions_;
		
	};
//---------------------------------------------------------------------------------------------------------------------

} // namespace GeometryProcess


#endif

