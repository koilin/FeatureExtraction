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

#include "MeshSegmentation.h"
#include <Core/PGMeshEntity.h>

#include <queue>
#include <time.h>
#include <fstream> 



namespace GeometryProcess
{
//---------------------------------------------------------------------------------------------------------------------

	// const std::string MeshSegmentation::STR_MESH_REGION = "__region__";
	const std::string MeshSegmentation::STR_SHARP_EDGE  = "__sharp_edge__";
	typedef MeshSegmentation::AdjacentGraph AdjacentGraph;
	typedef MeshSegmentation::SegmentationParameters SegmentationParameters;
	typedef MeshSegmentation::MeshRegions MeshRegions;

	MeshSegmentation::MeshSegmentation()
	{
		mesh_entity_ = NULL;
		n_regions_ = 0;
		param_.normal_variance_ = 0.9;
	}


	MeshSegmentation::~MeshSegmentation()
	{
	}

	void MeshSegmentation::set_mesh_entity(BaseMeshEntity * _mesh_entity)
	{
		mesh_entity_ = _mesh_entity;
	}
	void MeshSegmentation::set_segmentation_parameter(const SegmentationParameters & _param)
	{
		param_ = _param;
	}
	SegmentationParameters MeshSegmentation::get_segmentation_parameter() const
	{
		return param_;

	}
	int MeshSegmentation::n_regions() const
	{
		return n_regions_;
	}
	const AdjacentGraph & MeshSegmentation::get_adjacent_graph() const
	{
		return adjacent_graph_;
	}
	void MeshSegmentation::segment_mesh()
	{
		if (mesh_entity_->get_mesh_type() == BaseMeshEntity::POLYGON_MESH)
		{
			clock_t start_time = clock();


			PGMesh * mesh;
			mesh = &((PGMeshEntity *)mesh_entity_)->get_mesh();


			//--- add the properties ---//
			OpenMesh::FPropHandleT<int> region_property;
			if (!mesh->get_property_handle(region_property, BaseMeshEntity::STR_REGION))
			{
				mesh->add_property(region_property, BaseMeshEntity::STR_REGION);
			}
			OpenMesh::EPropHandleT<bool> sharp_edge_property;
			if (!mesh->get_property_handle(sharp_edge_property, STR_SHARP_EDGE))
			{
				mesh->add_property(sharp_edge_property, STR_SHARP_EDGE);
			}

			mesh->request_face_normals();
			mesh->request_vertex_normals();
			mesh->update_normals();

			//----------------find the sharp edge-------------------------//
			PGMesh::EdgeIter e_it;
			PGMesh::EdgeIter e_end(mesh->edges_end());
			double dihedral_angle;
			for (e_it = mesh->edges_begin(); e_it != e_end; ++e_it)
			{
				dihedral_angle = fabs(mesh->calc_dihedral_angle_fast(*e_it) * 180.0 / M_PI);
				if (dihedral_angle > param_.normal_variance_)
				{
					mesh->property(sharp_edge_property, *e_it) = true;
				}
				else
				{
					mesh->property(sharp_edge_property, *e_it) = false;
				}
			} // end of for


			unsigned int region_id;
			region_id = 0;
			PGMesh::FaceIter f_it;
			PGMesh::FaceIter f_end(mesh->faces_end());
			std::set<PGMesh::FaceHandle> faces;
			for (f_it = mesh->faces_begin(); f_it != f_end; ++f_it)
			{
				mesh->property(region_property, *f_it) = -1;
				faces.insert(*f_it);
			}

			//---------------------segment the mesh using bfs algorithm-------------//
			std::queue<PGMesh::FaceHandle> bfs;
			unsigned int region;
			region = 0;
			PGMesh::FaceHandle cur;
			std::set<PGMesh::FaceHandle>::iterator it;
			while (!faces.empty())
			{
				bfs.push(*faces.begin());
				mesh->property(region_property, *faces.begin()) = region;
				faces.erase(faces.begin());
				while (!bfs.empty())
				{
					cur = bfs.front();
					bfs.pop();
					for (PGMesh::FaceEdgeIter fe_it = mesh->fe_iter(cur); fe_it.is_valid(); ++fe_it)
					{
						if (!mesh->property(sharp_edge_property, *fe_it))
						{
							PGMesh::HalfedgeHandle & hh = mesh->halfedge_handle(*fe_it, 0);
							PGMesh::FaceHandle fh = mesh->face_handle(hh);
							if (fh == cur)
							{
								fh = mesh->face_handle(mesh->opposite_halfedge_handle(hh));
							}
							if (mesh->property(region_property, fh) == -1)
							{
								bfs.push(fh);
								mesh->property(region_property, fh) = region;
								it = faces.find(fh);
								if (it != faces.end())
								{
									faces.erase(it);
								}
							}
						}
					} // end of for 
				} // end of while bfs

				++ region;

			} // end of while faces
			n_regions_ = region + 1;

			//--- extract the adjacent graph of the regions ---//
			/*adjacent_graph_.clear();
			for (f_it = mesh->faces_begin(); f_it != f_end; ++f_it)
			{
				int &id = mesh->property(region_property, *f_it);
				for (PGMesh::FaceFaceIter ff_it = mesh->ff_iter(*f_it); ff_it.is_valid(); ++ ff_it)
				{
					int &aid = mesh->property(region_property, *ff_it);
					if (aid != id)
					{
						adjacent_graph_[id].insert(aid);
					}
				}
			} // end of for
			*/




			clock_t end_time = clock();

			std::ofstream f("D:\\codes\\GeometryProcessFramework\\time.txt");
			f << end_time - start_time << std::endl;

			f.close();



		} // end of if 

	}
	MeshRegions MeshSegmentation::segment_mesh_entity()
	{


		PGMesh * mesh;
		mesh = &((PGMeshEntity *)mesh_entity_)->get_mesh();

		clock_t start_time = clock();

		MeshRegions regions;


		//--- add the properties ---//
		OpenMesh::FPropHandleT<int> region_property;
		if (!mesh->get_property_handle(region_property, BaseMeshEntity::STR_REGION))
		{
			mesh->add_property(region_property, BaseMeshEntity::STR_REGION);
		}

		mesh->request_face_normals();
		mesh->update_face_normals();

		//------------ segment the mesh ---------------------------//
		int fn = mesh->n_faces();
		std::vector<bool> visited(fn, false);

		OpenMesh::Vec3d avgn;
		OpenMesh::Vec3d curn;
		OpenMesh::Vec3d totaln;
		double variance;

		int idx; // face index
		idx = 0;
		int rid; // region id
		rid = 0;
		int rn;  // region number
		rn = 0;

		std::queue<PGMesh::FaceHandle> Q;
		
		//++ idx;
		

		PGMesh::FaceHandle fh;

		while (idx < fn)
		{
			std::vector<PGMesh::FaceHandle> faces;

			for (; idx < fn; ++idx)
			{
				if (!visited[idx])
				{
					Q.push(PGMesh::FaceHandle(idx));
					faces.push_back(PGMesh::FaceHandle(idx));
					break;
				}
			}

			if (idx >= fn)
			{
				break;
			}

			mesh->property(region_property, PGMesh::FaceHandle(idx)) = rid;
			visited[idx] = true;

			avgn = mesh->normal(PGMesh::FaceHandle(idx));
			totaln = avgn;
			rn = 1;

			while (!Q.empty())
			{
				fh = Q.front();
				Q.pop();

				for (PGMesh::FaceFaceIter ff_it = mesh->ff_iter(fh); ff_it.is_valid(); ++ff_it)
				{
					int fid = (*ff_it).idx();
					if (!visited[(*ff_it).idx()])
					{
						variance = mesh->normal(*ff_it) | avgn;
						if (variance > param_.normal_variance_)
						{
							++ rn;
							Q.push(*ff_it);
							faces.push_back(*ff_it);
							mesh->property(region_property, *ff_it) = rid;
							visited[(*ff_it).idx()] = true;	

							totaln += (OpenMesh::Vec3d)mesh->normal(*ff_it);
							avgn = totaln / rn;
						}
					}
				}
			}
			if (!faces.empty())
			{
				regions.push_back(faces);
			}

			++ rid;
		}

		n_regions_ = rid;
		
		



		//---------------------------------------------------------//

		clock_t end_time = clock();

		std::ofstream f("D:\\codes\\GeometryProcessFramework\\time.txt");
		f << end_time - start_time << std::endl;

		f.close();

		return regions;

	}
	

//---------------------------------------------------------------------------------------------------------------------
} // namespace GeometryProcess
