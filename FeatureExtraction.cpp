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



#include "FeatureExtraction.h"
#include "MeshSegmentation.h"
#include "GeometryMath.h"
#include "LinearAlgebra.h"
#include "PCA.h"
#include <time.h>

#define PLANE_EPS 1e-6

using namespace GeometryProcess::Algorithm::LinearAgebra;
using namespace Math;

namespace GeometryProcess
{
//---------------------------------------------------------------------------------------------------------------------
	FeatureExtraction::FeatureExtraction()
	{
		mesh_entity_ = NULL;
	}


	FeatureExtraction::~FeatureExtraction()
	{
	}
	void FeatureExtraction::set_mesh_entity(BaseMeshEntity * _mesh_entity)
	{
		mesh_entity_ = _mesh_entity;
	}
	void FeatureExtraction::extract_features()
	{
		//OpenMesh::Vec3d top_left;
		//OpenMesh::Vec3d top_right;
		//OpenMesh::Vec3d bottom_left;
		//OpenMesh::Vec3d bottom_right;
		//extract_feature_points(top_left_, top_right_, bottom_right_, bottom_left_);

		
		if (!mesh_entity_)
		{
			return;
		}
		PGMesh * mesh;
		mesh = &((PGMeshEntity *)mesh_entity_)->get_mesh();


		clock_t start_time = clock();

		OpenMesh::Vec3d dirs[3];
		OpenMesh::Vec3d center;
		((PGMeshEntity *)mesh_entity_)->compute_PCA(center, dirs[0], dirs[1], dirs[2]);
		
		//--- segment the mesh ---//

		MeshSegmentation segment;
		segment.set_mesh_entity(mesh_entity_);
		MeshSegmentation::MeshRegions regions;
		regions = segment.segment_mesh_entity();

		//--- choose top n regions by area ---//

		std::set<int> indices;
		const int top_n = (10 > regions.size()) ? regions.size() : 10;

		double selected_idx;
		double k_max_region;

		for (int i = 0; i < top_n; ++i)
		{
			selected_idx = 0;
			k_max_region = 0;

			for (int k = 0; k < regions.size(); ++k)
			{
				if ((indices.find(k) == indices.end()) && (regions[k].region_area_ > k_max_region))
				{
					selected_idx = k;
					k_max_region = regions[k].region_area_;
				}
			}
			indices.insert(selected_idx);
			fiting_region_plane(regions[selected_idx]);
		}

		int bottom_region_idx;
		int top_region_idx;
		int inverse_mesh;
		bool swap_flag;
		swap_flag = false;
		inverse_mesh = extract_feature_points(center, dirs, regions, indices, bottom_region_idx, top_region_idx, swap_flag);

		//--- swap the second and the third directions ---//
		if (swap_flag)
		{
			std::swap(dirs[1], dirs[2]);
			inverse_mesh = extract_feature_points(center, dirs, regions, indices, bottom_region_idx, top_region_idx, swap_flag);
		}

		clock_t end_time = clock();

		recify_mesh_entity(center, dirs, inverse_mesh);
		

		//--- for debug, output four PCA points ---//

		std::ofstream f("running_time.txt");

		f << end_time - start_time << " ms" << std::endl;

		f.close();


		

	}
	//void FeatureExtraction::extract_feature_points(OpenMesh::Vec3d & _top_left, OpenMesh::Vec3d & _top_right, OpenMesh::Vec3d & _bottom_right, OpenMesh::Vec3d & _bottom_left)
	//{
	//	if (!mesh_entity_)
	//	{
	//		return;
	//	}

	//	if (mesh_entity_->get_mesh_type() != BaseMeshEntity::POLYGON_MESH)
	//	{
	//		return;
	//	}
	//	PGMesh * mesh;
	//	mesh = &((PGMeshEntity *)mesh_entity_)->get_mesh();

	//	clock_t start_time = clock();


	//	OpenMesh::Vec3d dirs[3];
	//	OpenMesh::Vec3d center;
	//	((PGMeshEntity *)mesh_entity_)->compute_PCA(center, dirs[0], dirs[1], dirs[2]);

	//	//--- if the mesh should be inversed ---//
	//	int inverse_mesh = 0;

	//	//--- segment the mesh ---//

	//	MeshSegmentation segment;
	//	segment.set_mesh_entity(mesh_entity_);
	//	//mesh->update_face_normals();
	//	MeshSegmentation::MeshRegions regions;
	//	regions = segment.segment_mesh_entity();

	//	//--- choose top n regions by area ---//

	//	std::set<int> indices;
	//	const int top_n = (10 > regions.size()) ? regions.size() : 10;

	//	double selected_idx;
	//	double k_max_region;

	//	for (int i = 0; i < top_n; ++i)
	//	{
	//		selected_idx = 0;
	//		k_max_region = 0;

	//		for (int k = 0; k < regions.size(); ++k)
	//		{
	//			if ((indices.find(k) == indices.end()) && (regions[k].region_area_ > k_max_region))
	//			{
	//				selected_idx = k;
	//				k_max_region = regions[k].region_area_;
	//			}
	//		}
	//		indices.insert(selected_idx);
	//		fiting_region_plane(regions[selected_idx]);
	//	}

	//	//--- find four PCA points ---//

	//	PGMesh::VertexHandle pca_left_ph;
	//	PGMesh::VertexHandle pca_right_ph;
	//	PGMesh::VertexHandle pca_top_ph;
	//	PGMesh::VertexHandle pca_bottom_ph;

	//	find_four_PCA_points(center, dirs, pca_left_ph, pca_right_ph, pca_top_ph, pca_bottom_ph);

	//	//--- find the bottom region ---//

	//	double max_region;
	//	max_region = 0;
	//	int bottom_idx;
	//	bottom_idx = 0;
	//	double cur;

	//	//OpenMesh::Vec3d y_axis(0, 1, 0);
	//	for (auto it = indices.begin(); it != indices.end(); ++it)
	//	{
	//		cur = fabs(regions[*it].plane_normal_ | dirs[1]) * regions[*it].region_area_;
	//		if (cur > max_region)
	//		{
	//			max_region = cur;
	//			bottom_idx = *it;
	//		}
	//	}

	//	//--- inverse the second direction ---//
	//	if ((regions[bottom_idx].avg_normal_ | dirs[1]) > 0)
	//	{
	//		std::swap(pca_top_ph, pca_bottom_ph);
	//		dirs[1] = -dirs[1];
	//		++inverse_mesh;
	//	}

	//	//--- compute the distance of the top point and the bottom region ---//
	//	const OpenMesh::Vec3d & tp = mesh->point(pca_top_ph);
	//	const OpenMesh::Vec3d & bp = mesh->point(pca_bottom_ph);
	//	const double last_height = fabs(tp[1] - bp[1]);
	//	double td;
	//	td = fabs((regions[bottom_idx].plane_normal_ | tp) + regions[bottom_idx].plane_d_);

	//	//--- the reference height, to filter some regions ---//
	//	const double ref_height = td / 2;

	//	//--- find the top region ---//
	//	int top_idx;
	//	top_idx = -1;

	//	double max_ref;
	//	max_ref = 0;
	//	double curdp;
	//	double cur_ref;

	//	auto bit = indices.find(bottom_idx);
	//	if (bit != indices.end())
	//	{
	//		//--- erase the index of the bottom region ---//
	//		indices.erase(bit);
	//	}
	//	for (auto it = indices.begin(); it != indices.end(); ++it)
	//	{
	//		//--- compute the distance of center to the fitting plane of the bottom region ---//
	//		curdp = fabs((regions[bottom_idx].plane_normal_ | regions[*it].center_) + regions[bottom_idx].plane_d_);
	//		if (curdp > ref_height)
	//		{
	//			//--- (normal variance) * (distance of the center to the fitting plane of the bottom region) ---//
	//			cur_ref = fabs(regions[bottom_idx].avg_normal_ | regions[*it].avg_normal_) * curdp;
	//			if (cur_ref > max_ref)
	//			{
	//				top_idx = *it;
	//				max_ref = cur_ref;
	//			}
	//		}
	//	}

	//	const OpenMesh::Vec3d & right_point = mesh->point(pca_right_ph);
	//	const OpenMesh::Vec3d & left_point = mesh->point(pca_left_ph);
	//	const double right_d = (right_point - regions[top_idx].center_).norm();
	//	const double left_d = (left_point - regions[top_idx].center_).norm();
	//	//--- swap the first direction ---//
	//	if (left_d > right_d)
	//	{
	//		std::swap(pca_left_ph, pca_right_ph);
	//		dirs[0] = -dirs[0];
	//		++inverse_mesh;
	//	}
	//	
	//	//--- get the bottom right feature point ---//
	//	ph_bottom_right_ = pca_right_ph;
	//	find_top_feature_points(regions[top_idx], ph_top_left_, ph_top_right_);
	//	find_bottom_left_point(regions[bottom_idx], ph_bottom_right_, ph_bottom_left_);

	//	clock_t end_time = clock();



	//	//-----------------------------------------------------------------------------------------//


	//	//--- for debug, output four PCA points ---//

	//	std::ofstream ffp("D:\\CWork\\codes\\GeometryProcessFramework\\four_feature_points.txt");

	//	ffp << ph_top_left_.idx() << std::endl;
	//	ffp << ph_top_right_.idx() << std::endl;
	//	ffp << ph_bottom_right_.idx() << std::endl;
	//	ffp << ph_bottom_left_.idx() << std::endl;

	//	ffp.close();




	//	

	//	//--- centering and reset the mesh to the PCA directions ---//
	//	OpenMesh::Vec3d cv;
	//	PGMesh::VertexIter v_end(mesh->vertices_end());
	//	double nx;
	//	double ny;
	//	double nz;
	//	for (auto v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
	//	{
	//		OpenMesh::Vec3d & p = mesh->point(*v_it);
	//		cv = p - center;
	//		nx = cv | dirs[0];
	//		ny = cv | dirs[1];
	//		nz = cv | dirs[2];
	//		mesh->set_point(*v_it, OpenMesh::Vec3d(nx, ny, nz));
	//	}

	//	mesh->update_normals();
	//	//if ((inverse_mesh % 2) == 1)
	//	//{
	//	//	((PGMeshEntity *)mesh_entity_)->inverse_mesh();
	//	//}
	//	((PGMeshEntity *)mesh_entity_)->update_bounding_box();


	//	//--- for debug, output four PCA points ---//

	//	std::ofstream fp("D:\\CWork\\codes\\GeometryProcessFramework\\four_key_points.txt");

	//	fp << pca_left_ph.idx() << std::endl;
	//	fp << pca_right_ph.idx() << std::endl;
	//	fp << pca_top_ph.idx() << std::endl;
	//	fp << pca_bottom_ph.idx() << std::endl;

	//	fp.close();


	//	//--- for debug ---//

	//	std::ofstream f("D:\\CWork\\codes\\GeometryProcessFramework\\max_selected.txt");

	//	f << "running time: " << end_time - start_time << std::endl;

	//	f << "height: " << last_height << std::endl;
	//	f << "top distance: " << td << std::endl;

	//	f << "Bottom region: \n";
	//	f << bottom_idx << std::endl;

	//	f << regions[bottom_idx].avg_normal_[0] << ", " << regions[bottom_idx].avg_normal_[1] << ", " << regions[bottom_idx].avg_normal_[2] << std::endl;

	//	f << "Top region: \n";
	//	f << top_idx << std::endl;
	//	f << regions[top_idx].region_area_ << std::endl;
	//	f << regions[top_idx].avg_normal_[0] << ", " << regions[top_idx].avg_normal_[1] << ", " << regions[top_idx].avg_normal_[2] << std::endl;
	//	f << regions[top_idx].plane_normal_[0] << ", " << regions[top_idx].plane_normal_[1] << ", " << regions[top_idx].plane_normal_[2] << std::endl;


	//	f.close();



	//}
	/** find the foure points after PCA adjustment */
	void FeatureExtraction::find_four_PCA_points(const OpenMesh::Vec3d & _c, const OpenMesh::Vec3d _dirs[3], PGMesh::VertexHandle & _left, PGMesh::VertexHandle & _right, PGMesh::VertexHandle & _top, PGMesh::VertexHandle & _bottom)
	{
		PGMesh * mesh = &((PGMeshEntity *)mesh_entity_)->get_mesh();
		double min_x;
		double max_x;
		min_x = 1e20;
		max_x = -1e20;

		double min_y;
		double max_y;
		min_y = 1e20;
		max_y = -1e20;

		double x;
		double y;

		PGMesh::VertexIter v_it;
		PGMesh::VertexIter v_end(mesh->vertices_end());

		for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
		{
			const OpenMesh::Vec3d & p = mesh->point(*v_it);
			x = (p - _c) | _dirs[0];
			y = (p - _c) | _dirs[1];

			if (x < min_x)
			{
				_left = *v_it;
				min_x = x;
			}
			if (x > max_x)
			{
				_right = *v_it;
				max_x = x;
			}
			if (y < min_y)
			{
				_bottom = *v_it;
				min_y = y;
			}
			if (y > max_y)
			{
				_top = *v_it;
				max_y = y;
			}
		}
		//std::ofstream f("D:\\CWork\\codes\\GeometryProcessFramework\\PCA_points.txt");

		//f << _left.idx() << std::endl;
		//f << _right.idx() << std::endl;
		//f << _bottom.idx() << std::endl;
		//f << _top.idx() << std::endl;


		//f.close();

	}

	/** fitting the plane of the input region */
	void FeatureExtraction::fiting_region_plane(MeshSegmentation::Region & _region)
	{
		std::set<int> phs;

		for (auto it = _region.faces_.begin(); it != _region.faces_.end(); ++it)
		{
			for (PGMesh::ConstFaceVertexIter cfv_it = _region.mesh_->cfv_begin(*it); cfv_it != _region.mesh_->cfv_end(*it); ++cfv_it)
			{
				phs.insert((*cfv_it).idx());
			}
		}

		std::vector<OpenMesh::Vec3d> points;
		for (auto pit = phs.begin(); pit != phs.end(); ++pit)
		{
			points.push_back(_region.mesh_->point(PGMesh::VertexHandle(*pit)));
		}
		fitting_plane(_region.plane_normal_, _region.plane_d_, points);
	}

	void FeatureExtraction::find_top_feature_points(const MeshSegmentation::Region & _region, PGMesh::VertexHandle & _top_left, PGMesh::VertexHandle & _top_right)
	{
		std::set<PGMesh::VertexHandle> pphs;
		PGMesh * mesh = _region.mesh_;

		for (auto it = _region.faces_.begin(); it != _region.faces_.end(); ++it)
		{
			for (auto cfv_it = mesh->cfv_iter(*it); cfv_it.is_valid(); ++cfv_it)
			{
				pphs.insert(*cfv_it);
			}
		}

		double d;

		PGMesh::ConstVertexIter v_end(mesh->vertices_end());

		for (auto v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
		{
			if (pphs.find(*v_it) == pphs.end())
			{
				const OpenMesh::Vec3d & p = mesh->point(*v_it);
				d = fabs((_region.plane_normal_ | p) + _region.plane_d_);
				if (d < PLANE_EPS)
				{
					pphs.insert(*v_it);
				}
			}
		}

		OpenMesh::Vec3d c;
		OpenMesh::Vec3d dirs[3];
		double * x;
		double * y;
		double * z;
		int nv;
		nv = pphs.size();
		x = new double[nv];
		y = new double[nv];
		z = new double[nv];

		int idx;
		idx = 0;

		for (auto vit = pphs.begin(); vit != pphs.end(); ++vit)
		{
			const OpenMesh::Vec3d & pp = mesh->point(*vit);
			x[idx] = pp[0];
			y[idx] = pp[1];
			z[idx] = pp[2];
			++idx;
		}
		PCA pca;
		pca.Compute(x, y, z, nv, &c[0], &dirs[0][0], &dirs[1][0], &dirs[2][0]);

		delete[] x;
		x = NULL;
		delete[] y;
		y = NULL;
		delete[] z;
		z = NULL;

		double min_x;
		double max_x;
		min_x = 1e20;
		max_x = -1e20;

		double px;
		for (auto it = pphs.begin(); it != pphs.end(); ++it)
		{
			const OpenMesh::Vec3d & p = mesh->point(*it);
			px = (p - c) | dirs[0];
			if (px < min_x)
			{
				min_x = px;
				_top_left = *it;
			}
			if (px > max_x)
			{
				max_x = px;
				_top_right = *it;
			}
		}
	}
	/** find the bottom left point */
	void FeatureExtraction::find_bottom_left_point(const MeshSegmentation::Region & _region, const PGMesh::VertexHandle & _bottom_right, PGMesh::VertexHandle & _bottom_left)
	{
		PGMesh * mesh = _region.mesh_;
		double max_ref;
		max_ref = -1e20;
		std::set<PGMesh::VertexHandle> pphs;

		for (auto it = _region.faces_.begin(); it != _region.faces_.end(); ++it)
		{
			for (auto cfv_it = mesh->cfv_iter(*it); cfv_it.is_valid(); ++cfv_it)
			{
				pphs.insert(*cfv_it);
			}
		}

		double d;
		const OpenMesh::Vec3d & brp = mesh->point(_bottom_right);
		for (auto vit = pphs.begin(); vit != pphs.end(); ++vit)
		{
			const OpenMesh::Vec3d & p = mesh->point(*vit);
			d = (brp - p).norm();
			if (d > max_ref)
			{
				max_ref = d;
				_bottom_left = *vit;
			}			
		}		
	}
	/** recify the mesh entity */
	void FeatureExtraction::recify_mesh_entity(const OpenMesh::Vec3d & _c, const OpenMesh::Vec3d _dirs[3], int _inverse_mesh)
	{
		PGMesh * mesh;
		mesh = &((PGMeshEntity *)mesh_entity_)->get_mesh();


		//--- centering and reset the mesh to the PCA directions ---//
		OpenMesh::Vec3d cv;
		PGMesh::VertexIter v_end(mesh->vertices_end());
		double nx;
		double ny;
		double nz;
		for (auto v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
		{
			OpenMesh::Vec3d & p = mesh->point(*v_it);
			cv = p - _c;
			nx = cv | _dirs[0];
			ny = cv | _dirs[1];
			nz = cv | _dirs[2];
			mesh->set_point(*v_it, OpenMesh::Vec3d(nx, ny, nz));
		}

		mesh->update_normals();
		//if ((inverse_mesh % 2) == 1)
		//{
		//	((PGMeshEntity *)mesh_entity_)->inverse_mesh();
		//}
		((PGMeshEntity *)mesh_entity_)->update_bounding_box();

	}

	int FeatureExtraction::extract_feature_points(const OpenMesh::Vec3d & _c, 
		                                          OpenMesh::Vec3d _dirs[3], 
		                                          const MeshSegmentation::MeshRegions & _regions, 
		                                          const std::set<int> & _top_n_regions,
		                                          int & _top_region_idx,
		                                          int & _bottom_region_idx, 
		                                          bool & _is_swap)
	{
		int inverse_mesh;
		inverse_mesh = 0;
		PGMesh * mesh;
		mesh = &((PGMeshEntity *)mesh_entity_)->get_mesh();


		//--- find four PCA points ---//

		PGMesh::VertexHandle pca_left_ph;
		PGMesh::VertexHandle pca_right_ph;
		PGMesh::VertexHandle pca_top_ph;
		PGMesh::VertexHandle pca_bottom_ph;

		find_four_PCA_points(_c, _dirs, pca_left_ph, pca_right_ph, pca_top_ph, pca_bottom_ph);
		//--- find the bottom region ---//

		double max_region;
		max_region = 0;
		int bottom_idx;
		bottom_idx = 0;
		double cur;

		//OpenMesh::Vec3d y_axis(0, 1, 0);
		for (auto it = _top_n_regions.begin(); it != _top_n_regions.end(); ++it)
		{
			cur = fabs(_regions[*it].plane_normal_ | _dirs[1]) * _regions[*it].region_area_;
			if (cur > max_region)
			{
				max_region = cur;
				bottom_idx = *it;
			}
		}

		//--- inverse the second direction ---//
		if ((_regions[bottom_idx].avg_normal_ | _dirs[1]) > 0)
		{
			std::swap(pca_top_ph, pca_bottom_ph);
			_dirs[1] = -(_dirs[1]);
			++ inverse_mesh;
		}
		//--- compute the distance of the top point and the bottom region ---//
		const OpenMesh::Vec3d & pca_tp = mesh->point(pca_top_ph);
		const OpenMesh::Vec3d & pca_bp = mesh->point(pca_bottom_ph);
		const double last_height =  fabs(((pca_tp - _c) | _dirs[1]) - ((pca_bp - _c) | _dirs[1]));
		double td;
		td = fabs((_regions[bottom_idx].plane_normal_ | pca_tp) + _regions[bottom_idx].plane_d_);

		//--- the reference height, to filter some regions ---//
		const double ref_height = td / 2;

		//--- find the top region ---//
		int top_idx;
		top_idx = -1;

		double max_ref;
		max_ref = 0;
		double curdp;
		double cur_ref;

		for (auto it = _top_n_regions.begin(); it != _top_n_regions.end(); ++it)
		{
			if ((*it) != bottom_idx)
			{
				//--- compute the distance of center to the fitting plane of the bottom region ---//
				curdp = fabs((_regions[bottom_idx].plane_normal_ | _regions[*it].center_) + _regions[bottom_idx].plane_d_);
				if (curdp > ref_height)
				{
					//--- (normal variance) * (distance of the center to the fitting plane of the bottom region) ---//
					cur_ref = fabs(_regions[bottom_idx].avg_normal_ | _regions[*it].avg_normal_) * curdp;
					if (cur_ref > max_ref)
					{
						top_idx = *it;
						max_ref = cur_ref;
					}
				}

			}
		} // end of for 

		const OpenMesh::Vec3d & right_point = mesh->point(pca_right_ph);
		const OpenMesh::Vec3d & left_point = mesh->point(pca_left_ph);
		const double right_d = (right_point - _regions[top_idx].center_).norm();
		const double left_d = (left_point - _regions[top_idx].center_).norm();
		//--- swap the first direction ---//
		if (left_d > right_d)
		{
			std::swap(pca_left_ph, pca_right_ph);
			_dirs[0] = -(_dirs[0]);
			++inverse_mesh;
		}

		//--- get the bottom right feature point ---//
		ph_bottom_right_ = pca_right_ph;
		find_top_feature_points(_regions[top_idx], ph_top_left_, ph_top_right_);
		find_bottom_left_point(_regions[bottom_idx], ph_bottom_right_, ph_bottom_left_);

		const OpenMesh::Vec3d & trp = mesh->point(ph_top_left_);
		const OpenMesh::Vec3d & brp = mesh->point(ph_bottom_left_);
		double d = fabs(((trp - _c) | _dirs[1]) - ((brp - _c) | _dirs[1]));


		if (d < (last_height * 0.5))
		{
			_is_swap = true;
		}





		//-----------------------------------------------------------------------------------------//


		//--- for debug, output four PCA points ---//

		std::ofstream ffp("four_feature_points.txt");

		ffp << ph_top_left_.idx() << std::endl;
		ffp << ph_top_right_.idx() << std::endl;
		ffp << ph_bottom_right_.idx() << std::endl;
		ffp << ph_bottom_left_.idx() << std::endl;

		ffp.close();




		/*

		//--- for debug, output four PCA points ---//

		std::ofstream fp("D:\\CWork\\codes\\GeometryProcessFramework\\four_key_points.txt");

		fp << pca_left_ph.idx() << std::endl;
		fp << pca_right_ph.idx() << std::endl;
		fp << pca_top_ph.idx() << std::endl;
		fp << pca_bottom_ph.idx() << std::endl;

		fp.close();


		//--- for debug ---//

		std::ofstream f("D:\\CWork\\codes\\GeometryProcessFramework\\max_selected.txt");

		f << "running time: " << end_time - start_time << std::endl;

		f << "height: " << last_height << std::endl;
		f << "top distance: " << td << std::endl;

		f << "Bottom region: \n";
		f << bottom_idx << std::endl;

		f << _regions[bottom_idx].avg_normal_[0] << ", " << _regions[bottom_idx].avg_normal_[1] << ", " << _regions[bottom_idx].avg_normal_[2] << std::endl;

		f << "Top region: \n";
		f << top_idx << std::endl;
		f << _regions[top_idx].region_area_ << std::endl;
		f << _regions[top_idx].avg_normal_[0] << ", " << _regions[top_idx].avg_normal_[1] << ", " << _regions[top_idx].avg_normal_[2] << std::endl;
		f << _regions[top_idx].plane_normal_[0] << ", " << _regions[top_idx].plane_normal_[1] << ", " << _regions[top_idx].plane_normal_[2] << std::endl;


		f.close();
		*/

		//----------------------------------------------------------------------------------------------------
		
		return inverse_mesh;
	}

//---------------------------------------------------------------------------------------------------------------------
} // namespace GeometryProcess
