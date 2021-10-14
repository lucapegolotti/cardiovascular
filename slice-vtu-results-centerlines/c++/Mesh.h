#include <iostream>
#include <set>
#include <cstring>
#include <string>
#include <vector>
#include <functional>
#include <map>

#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

class Graphics;

#ifndef MESH_H
#define MESH_H

class Mesh {
  public:
    Mesh();
    ~Mesh();

    void add_geometry();
    void compute_plane_dist(double position[3], double normal[3]);
    void extract_all_slices(vtkPolyData* centerlines,
                            bool compute_average_fields = false,
                            bool update_graphics = true);
    void extract_slice(double pos[3], double inscribedRadius, double normal[3]);
    vtkSmartPointer<vtkPolyData> find_best_slice(double position[3], vtkPolyData* isosurface);
    void read_mesh(const std::string& fileName);
    void remove_data_arrays(const std::set<std::string>& slice_data_names);
    vtkSmartPointer<vtkPolyData> trim_slice(vtkPolyData* slice, double position[3], double radius);
    void write_vtp(vtkPolyData* slice, std::string filename);
    void write_centerlines_and_fields(vtkPolyData* slice, std::string filename);
    double integrate_on_slice(vtkPolyData* slice, vtkIdType numcells,
                              std::vector<double> area_cells,
                              std::function<double(vtkIdType)> fun);
    double compute_area_slice(vtkPolyData* slice, std::vector<double>& area_cells);
    Graphics* graphics_;

    std::string mesh_file_name_;

    // The name of the scalar field used to slice the mesh.
    const char* slice_scalar_name_;

    vtkSmartPointer<vtkDoubleArray> plane_dist_;

    vtkSmartPointer<vtkPolyData> mesh_polydata_;

    bool trim_slice_using_incribed_sphere_;

    vtkSmartPointer<vtkUnstructuredGrid> unstructured_mesh_;

    // avg pressures over each slice (key: name in the original vtu)
    std::map<std::string,vtkSmartPointer<vtkDoubleArray>> avg_pressures_;

    // flowrate over each slice (key: name in the original vtu)
    std::map<std::string,vtkSmartPointer<vtkDoubleArray>> flowrates_;

    vtkSmartPointer<vtkDoubleArray> areas_;
};

#endif
