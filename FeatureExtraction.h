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


#ifndef _FEATURE_EXTRACTION_H_
#define _FEATURE_EXTRACTION_H_

#include <map>
#include <set>
#include <vector>
#include <Core/BaseMeshEntity.h>
#include <Core/PGMeshEntity.h>
#include "MeshSegmentation.h"



namespace GeometryProcess
{
	class FeatureExtraction
	{
	public:
		typedef std::vector<OpenMesh::Vec3d> FeatureLine;
	public:
		FeatureExtraction();
		~FeatureExtraction();
	public:
		void set_mesh_entity(BaseMeshEntity * _mesh_entity);
		void extract_features();
		int extract_feature_points(const OpenMesh::Vec3d & _c, 
			                       OpenMesh::Vec3d _dirs[3], 
			                       const MeshSegmentation::MeshRegions & _regions, 
			                       const std::set<int> & _top_n_regions,
			                       int & _top_region_idx,
			                       int & _bottom_region_idx, 
			                       bool & _swap);
		/** generate the recified directions for the mesh entity */
		void generate_recify_directions(OpenMesh::Vec3d & _center, OpenMesh::Vec3d _recified_dirs[3], const OpenMesh::Vec3d _pca_dirs[3], const MeshSegmentation::Region & _bottom_region);
	private:
		/** find the foure points after PCA adjustment */
		void find_four_PCA_points(const OpenMesh::Vec3d & _c, const OpenMesh::Vec3d _dirs[3], PGMesh::VertexHandle & _left, PGMesh::VertexHandle & _right, PGMesh::VertexHandle & _top, PGMesh::VertexHandle & _bottom);
		/** fitting the plane of the input region */
		void fiting_region_plane(MeshSegmentation::Region & _region);
		/** find two features of the top region */
		void find_top_feature_points(const MeshSegmentation::Region & _region, PGMesh::VertexHandle & _top_left, PGMesh::VertexHandle & _top_right);
		/** find the bottom left point */
		void find_bottom_left_point(const MeshSegmentation::Region & _region, const PGMesh::VertexHandle & _bottom_right, PGMesh::VertexHandle & _bottom_left);
		/** recify the mesh entity */
		void recify_mesh_entity(const OpenMesh::Vec3d & _c, const OpenMesh::Vec3d _dirs[3], int _inverse_mesh = 0);
	private:
		/** extract the left feature line */
		void extract_left_line();
	private:
		BaseMeshEntity * mesh_entity_;
		PGMesh::VertexHandle ph_top_left_;
		PGMesh::VertexHandle ph_top_right_;
		PGMesh::VertexHandle ph_bottom_left_;
		PGMesh::VertexHandle ph_bottom_right_;
		OpenMesh::Vec3d top_left_point_;
		OpenMesh::Vec3d top_right_point_;
		OpenMesh::Vec3d bottom_left_point_;
		OpenMesh::Vec3d bottom_right_point_;
		FeatureLine left_line_;
	};
}


#endif // !_FEATURE_EXTRACTION_H_
