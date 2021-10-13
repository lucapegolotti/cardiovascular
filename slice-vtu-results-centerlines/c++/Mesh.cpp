
#include "Slice.h"
#include "Mesh.h"

#include <vtkCleanPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkContourGrid.h>
#include <vtkDoubleArray.h>
#include <vtkExtractGeometry.h>
#include <vtkGeometryFilter.h>
#include <vtkGenericCell.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkSphere.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkPolyDataWriter.h>

#include <iostream>
#include <string>
#include <chrono>
#include <math.h>

//------
// Mesh
//------
//
Mesh::Mesh()
{
  slice_scalar_name_ = "plane_dist";

  trim_slice_using_incribed_sphere_ = false;
}

Mesh::~Mesh()
{
}

//-----------
// read_mesh
//-----------
// Read a volume mesh from a vtu file.
//
void Mesh::read_mesh(const std::string& file_name)
{
  mesh_file_name_ = file_name;
  std::cout << "[read_mesh] Read mesh: " << file_name << std::endl;
  unstructured_mesh_ = vtkSmartPointer<vtkUnstructuredGrid>::New();

  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();
  unstructured_mesh_->DeepCopy(reader->GetOutput());
  std::cout << "[read_mesh] Number of elements: " << unstructured_mesh_->GetNumberOfCells() << std::endl;

  // Remove data arrays we don't want to generate slice data for.
  std::set<std::string> retain_data_names = { "pressure", "velocity" };   // Data to not remove.
  remove_data_arrays(retain_data_names);

  // Convert the unstructured mesh to polydata for visualization.
  auto geometry_filter = vtkSmartPointer<vtkGeometryFilter>::New();
  geometry_filter->SetInputData(reader->GetOutput());
  geometry_filter->Update();
  mesh_polydata_ = vtkSmartPointer<vtkPolyData>::New();
  mesh_polydata_->DeepCopy(geometry_filter->GetOutput());

  // Add a point data array to store plane distance.
  int num_pts = unstructured_mesh_->GetNumberOfPoints();
  plane_dist_ = vtkSmartPointer<vtkDoubleArray>::New();
  plane_dist_->SetName(slice_scalar_name_);
  plane_dist_->SetNumberOfComponents(1);
  plane_dist_->SetNumberOfTuples(num_pts);
  for (int i = 0; i < num_pts; i++) {
    plane_dist_->SetValue(i,0.0);
  }
  unstructured_mesh_->GetPointData()->AddArray(plane_dist_);

  std::cout << "[read_mesh] Done. " << std::endl;
}

//--------------
// add_geometry
//--------------
// Add the mesh geometry for visualization.
//
void Mesh::add_geometry()
{
  auto geom = graphics_->create_geometry(mesh_polydata_);
  geom->GetProperty()->SetColor(0.8, 0.0, 0.0);
  //geom->GetProperty()->SetRepresentationToWireframe();
  geom->GetProperty()->EdgeVisibilityOn();
  geom->GetProperty()->SetColor(0.8, 0.8, 0.8);
  geom->GetProperty()->SetOpacity(0.3);
  //geom->GetProperty()->SetRepresentationToWireframe();
  //geom->GetProperty()->FrontfaceCullingOn();
  //geom->GetProperty()->BackfaceCullingOn();
  geom->PickableOff();
  graphics_->add_geometry(geom);
}

//--------------------
// integrate_on_slice
//--------------------
// Integrate function fun on a slice. Function fun must be evaluated at the index of
// a slice point
//
double Mesh::integrate_on_slice(vtkPolyData* slice, vtkIdType numcells,
                                double area_cells[],
                                std::function<double(vtkIdType)> fun)
{
  double sum = 0;
  for (int j = 0; j < numcells; j++)
  {
    vtkIdType cell_size;
    const vtkIdType* cell_points;
    slice->GetPolys()->GetCellAtId(j, cell_size, cell_points);
    double cur_sum = 0;
    for (int k = 0; k < cell_size; k++)
      cur_sum  += fun(cell_points[k]);
    sum += area_cells[j] * (cur_sum) / 3.;
  }
  return sum;
}

//--------------------
// extract_all_slices
//--------------------
// Extract slices for all centerline points.
//
void Mesh::extract_all_slices(vtkPolyData* centerlines, bool compute_average_fields, bool update_graphics)
{
  std::cout << "[extract_all_slices] " << std::endl;
  std::cout << "[extract_all_slices] ========== Mesh::extract_all_slices ========== " << std::endl;
  auto start_time = std::chrono::steady_clock::now();

  // Set the scalar field used to extract a slice.
  unstructured_mesh_->GetPointData()->SetActiveScalars(slice_scalar_name_);

  // Get centerline data.
  int num_points = centerlines->GetNumberOfPoints();
  auto points = centerlines->GetPoints();
  auto normal_data = vtkDoubleArray::SafeDownCast(centerlines->GetPointData()->GetArray("CenterlineSectionNormal"));
  auto radius_data = vtkDoubleArray::SafeDownCast(centerlines->GetPointData()->GetArray("MaximumInscribedSphereRadius"));
  std::cout << "[extract_all_slices] Number of centerline points: " << num_points << std::endl;

  if (compute_average_fields) {
      areas_ = vtkSmartPointer<vtkDoubleArray>::New();
      areas_->SetName("area");
      areas_->SetNumberOfComponents(1);
      areas_->SetNumberOfTuples(num_points);
  }

  // Extract slices.
  //
  for (int i = 0; i < num_points; i++) {
    std::cout << "i/num_points = " << i << "/" << num_points << std::endl << std::flush;

    double position[3];
    double normal[3];
    points->GetPoint(i, position);
    normal_data->GetTuple(i, normal);
    double radius = radius_data->GetValue(i);
    // Compute the distance of each mesh node to the plane.
    compute_plane_dist(position, normal);

    // Extract the slice.
    auto contour = vtkSmartPointer<vtkContourGrid>::New();
    contour->SetInputData(unstructured_mesh_);
    contour->SetValue(0, 0.0);
    contour->ComputeNormalsOff();
    contour->Update();
    vtkSmartPointer<vtkPolyData> slice = vtkSmartPointer<vtkPolyData>::New();
    slice->DeepCopy(contour->GetOutput());
    // Trim the slice using the incribed sphere radius.
    std::cout << 1 << std::endl << std::flush;
    if (trim_slice_using_incribed_sphere_) {
      slice = trim_slice(slice.GetPointer(), position, radius);
      std::cout << 2 << std::endl << std::flush;
    // Find the slice geometry center closest to the selected point.
    } else {
      slice = find_best_slice(position, slice.GetPointer());
    }
    if (compute_average_fields) {
      vtkIdType num_point_arrays = slice->GetPointData()->GetNumberOfArrays();
      int num_points_slice = slice->GetNumberOfPoints();
      auto polys = slice->GetPolys();
      int numcells = polys->GetNumberOfCells();
      auto points_slice = slice->GetPoints();
      // std::stringstream file_name;
      // file_name << "slice-1.vtp";
      // write_vtp(slice, file_name.str());

      double area_cells[numcells];
      double area = 0;
      for (int icell = 0; icell < numcells; icell++)
      {
  vtkIdType cell_size;
  const vtkIdType* cell_points;
  polys->GetCellAtId(icell, cell_size, cell_points);
  double edge_sizes[cell_size];

  // from now on, we assume that the elements are triangles (true for
  // tetrahedral mesh)

  double nodes[cell_size][3];
  for (int j = 0; j < cell_size; j++)
  {
    double* points = slice->GetPoints()->GetPoint(cell_points[j]);
    for (int k = 0; k < 3; k++)
      nodes[j][k] = points[k];
  }
  for (int j = 0; j < cell_size; j++)
  {
    double diff[3];
    for (int k = 0; k < 3; k++)
      diff[k] = nodes[(j + 1) % 3][k] - nodes[j][k];

    double sum = 0;
    for (int k = 0; k < 3; k++)
      sum += diff[k] * diff[k];
    edge_sizes[j] = std::sqrt(sum);
  }

  double s = (edge_sizes[0] + edge_sizes[1] + edge_sizes[2]) / 2;
  area_cells[icell] = std::sqrt(s * (s - edge_sizes[0]) *
             (s - edge_sizes[1]) *
             (s - edge_sizes[2]));
  area += area_cells[icell];
      }
  areas_->SetValue(i, area);
  for (int ipoint = 0; ipoint < num_point_arrays; ipoint++) {
  int type = slice->GetPointData()->GetArray(ipoint)->GetDataType();
  auto name = std::string(slice->GetPointData()->GetArrayName(ipoint));
  int upos = name.find('_');
  if (std::strcmp(name.substr(0, upos).c_str(),"velocity") == 0) {
    double flux[num_points_slice];
    for (int j = 0; j < num_points_slice; j++) {
      double* vel = slice->GetPointData()->GetArray(ipoint)->GetTuple3(j);
      double curflux = 0;
      for (int k = 0; k < 3; k++) {
        curflux += vel[k] * normal[k];
      }
      flux[j] = curflux;
    }
    double integral = integrate_on_slice(slice, numcells, area_cells,
            [&](vtkIdType index) -> double {
              return flux[index];
            });
    if (flowrates_.find(name) == flowrates_.end()) {
      flowrates_[name] = vtkSmartPointer<vtkDoubleArray>::New();
      flowrates_[name]->SetName(name.c_str());
      flowrates_[name]->SetNumberOfComponents(1);
      flowrates_[name]->SetNumberOfTuples(num_points);
    }
    flowrates_[name]->SetValue(i, integral/area);
  }
  if (std::strcmp(name.substr(0, upos).c_str(),"pressure") == 0) {
    double integral = integrate_on_slice(slice, numcells, area_cells,
            [&](vtkIdType index) -> double {
              return slice->GetPointData()->GetArray(ipoint)->GetTuple1(index);
            });
    if (avg_pressures_.find(name) == avg_pressures_.end()) {
      avg_pressures_[name] = vtkSmartPointer<vtkDoubleArray>::New();
      avg_pressures_[name]->SetName(name.c_str());
      avg_pressures_[name]->SetNumberOfComponents(1);
      avg_pressures_[name]->SetNumberOfTuples(num_points);
    }
    avg_pressures_[name]->SetValue(i, integral/area);
  }
      }
      // slice->Delete();
    }

    if (update_graphics) {
      auto geom = graphics_->create_geometry(slice);
      geom->GetProperty()->SetColor(0.8, 0.0, 0.0);
      geom->PickableOff();
      graphics_->add_geometry(geom);
    }
  }

  auto end_time = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end_time - start_time;
  std::cout << "[extract_all_slices] Elapsed time: " << elapsed_seconds.count() << "s\n";
}

//--------------------
// remove_data_arrays
//--------------------
// Remove data arrays we don't want to generate slice data for.
//
// Arguments:
//   retain_data_names - The data names to not to be removed.
//
void Mesh::remove_data_arrays(const std::set<std::string>& retain_data_names)
{
  std::cout << "[remove_data_arrays] " << std::endl;
  std::cout << "[remove_data_arrays] ========== Mesh::remove_data_arrays ========== " << std::endl;

  // Note: Can't remove arrays in this loop because it crashes vtk,
  // changes unstructured_mesh_ data it seems.
  //
  vtkIdType num_point_arrays = unstructured_mesh_->GetPointData()->GetNumberOfArrays();
  std::cout << "[remove_data_arrays] Number of node data arrays: " << num_point_arrays << std::endl;
  std::vector<std::string> remove_data_names;

  for (int i = 0; i < num_point_arrays; i++) {
    int type = unstructured_mesh_->GetPointData()->GetArray(i)->GetDataType();
    auto name = std::string(unstructured_mesh_->GetPointData()->GetArrayName(i));
    int upos = name.find('_');
    int num_comp = unstructured_mesh_->GetPointData()->GetNumberOfComponents();
    if (retain_data_names.count(name.substr(0, upos).c_str()) == 0) {
      std::cout << "[remove_data_arrays] --- Remove data array: " << name << std::endl;
      remove_data_names.push_back(name);
    } else {
      std::cout << "[remove_data_arrays] +++ Retain data array: " << name << std::endl;
    }
  }

  for (auto const& name : remove_data_names) {
    unstructured_mesh_->GetPointData()->RemoveArray(name.c_str());
  }
}

//--------------------
// compute_plane_dist
//--------------------
// Compute the distance from each point in the mesh to a plane.
//
// The plane distance is store in the plane_dist_ Point Data array.
//
// Arguments:
//   position - The plane position.
//   normal - The plane normal.
//
void Mesh::compute_plane_dist(double position[3], double normal[3])
{
  int num_pts = unstructured_mesh_->GetNumberOfPoints();
  auto mesh_points = unstructured_mesh_->GetPoints();

  // Compute distance of each mesh point from the slicing plane.
  for (int i = 0; i < num_pts; i++) {
    double pt[3];
    mesh_points->GetPoint(i, pt);
    double d = normal[0] * (pt[0] - position[0]) +
         normal[1] * (pt[1] - position[1]) +
         normal[2] * (pt[2] - position[2]);
    plane_dist_->SetValue(i, d);
  }
}

//----------------
// extract_slice
//----------------
// Extract a slice oriented with the given normal at the given position.
//
void Mesh::extract_slice(double position[3], double inscribedRadius, double normal[3])
{
  std::cout << "[extract_slice] " << std::endl;
  std::cout << "[extract_slice] ========== Mesh::extract_slice ========== " << std::endl;
  int num_pts = unstructured_mesh_->GetNumberOfPoints();
  std::cout << "[extract_slice] num_pts: " << num_pts << std::endl;
  double pt[3];
  auto mesh_points = unstructured_mesh_->GetPoints();
  std::cout << "[extract_slice] Normal: " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;

  // Compute distance of each mesh point from the slicing plane.
  compute_plane_dist(position, normal);

  // Extract isosurface.
  unstructured_mesh_->GetPointData()->SetActiveScalars(slice_scalar_name_);
  auto contour = vtkSmartPointer<vtkContourGrid>::New();
  contour->SetInputData(unstructured_mesh_);
  contour->SetValue(0, 0.0);
  contour->ComputeScalarsOn();
  contour->Update();
  auto slice = contour->GetOutput();
  int num_slice_pts = slice ->GetNumberOfPoints();
  std::cout << "[extract_slice] num_slice_pts: " << num_slice_pts << std::endl;

  // Show the slice.
  /*
  auto slice_geom = graphics_->create_geometry(slice);
  slice_geom->GetProperty()->SetRepresentationToWireframe();
  slice_geom->GetProperty()->SetColor(0.0, 0.5, 0.0);
  slice_geom->GetProperty()->SetLineWidth(1);
  slice_geom->PickableOff();
  graphics_->add_geometry(slice_geom);
  */

  // Trim the slice using the incribed sphere radius.
  if (trim_slice_using_incribed_sphere_) {
    slice = trim_slice(slice, position, inscribedRadius);

  // For multipled disjoint slice geometry find the geometry with
  // center closest to the selected point.
  } else {
    slice = find_best_slice(position, slice);
  }

  // Write the slice to a .vtp file.
  std::stringstream file_name;
  file_name << "slice-1.vtp";
  write_vtp(slice, file_name.str());

  // Show the trimmed slice.
  auto geom = graphics_->create_geometry(slice);
  geom->GetProperty()->SetColor(1.0, 0.0, 0.0);
  geom->GetProperty()->EdgeVisibilityOn();
  geom->PickableOff();
  graphics_->add_geometry(geom);
}

//------------
// trim_slice
//------------
// Trim the slice using the incribed sphere radius.
//
vtkSmartPointer<vtkPolyData>
Mesh::trim_slice(vtkPolyData* slice, double position[3], double radius)
{
  auto points = slice->GetPoints();
  int num_points = slice->GetNumberOfPoints();
  int num_cells = slice->GetNumberOfCells();
  double center[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < num_points; i++) {
    double pt[3];
    points->GetPoint(i, pt);
    center[0] += pt[0];
    center[1] += pt[1];
    center[2] += pt[2];
  }

  center[0] /= num_points;
  center[1] /= num_points;
  center[2] /= num_points;

  // Creae a sphere implicit function.
  auto sphere = vtkSmartPointer<vtkSphere>::New();
  sphere->SetCenter(position);
  sphere->SetRadius(radius);
  // Extract geometry within the sphere implicit function.
  //
  std::cout << ".1" << std::endl << std::flush;
  vtkSmartPointer<vtkPolyData> trimmed_slice;
 
  if (false) {
    auto extract = vtkSmartPointer<vtkExtractGeometry>::New();
    extract->SetInputData(slice);
    extract->SetImplicitFunction(sphere);
    extract->ExtractInsideOn();
    extract->ExtractBoundaryCellsOn();
    extract->Update();
    auto extract_grid = extract->GetOutput();

    std::cout << ".2" << std::endl << std::flush;
    auto geometry_filter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometry_filter->SetInputData(extract_grid);
    geometry_filter->Update();

    trimmed_slice = vtkSmartPointer<vtkPolyData>::New();
    trimmed_slice->DeepCopy(geometry_filter->GetOutput());

  // Clip geometry to the sphere implicit function.

  } else {
    auto clip = vtkSmartPointer<vtkClipPolyData>::New();
    clip->SetInputData(slice);
    clip->SetClipFunction(sphere);
    clip->InsideOutOn();
    clip->Update();
    std::cout << "..2" << std::endl << std::flush;
    trimmed_slice = vtkSmartPointer<vtkPolyData>::New();
    trimmed_slice->DeepCopy(clip->GetOutput());
    std::cout << "..3" << std::endl << std::flush;
    }
  
  return trimmed_slice;
}

//-----------------
// find_best_slice
//-----------------
// For multipled disjoint slice geometry find the geometry with center
// closest to the selected point.
//
vtkSmartPointer<vtkPolyData>
Mesh::find_best_slice(double position[3], vtkPolyData* slice)
{
  auto conn_filter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  conn_filter->SetInputData(slice);
  conn_filter->SetExtractionModeToSpecifiedRegions();

  double min_d = 1e6;
  int rid = 0;
  double center[3] = {0.0, 0.0, 0.0};
  vtkSmartPointer<vtkPolyData> min_comp;

  while (true) {
    conn_filter->AddSpecifiedRegion(rid);
    conn_filter->Update();
    auto component = vtkSmartPointer<vtkPolyData>::New();
    component->DeepCopy(conn_filter->GetOutput());
    if (component->GetNumberOfCells() <= 0) {
      break;
    }
    conn_filter->DeleteSpecifiedRegion(rid);

    auto clean_filter = vtkSmartPointer<vtkCleanPolyData>::New();
    clean_filter->SetInputData(component);
    clean_filter->Update();
    component->DeepCopy(clean_filter->GetOutput());

    auto comp_points = component->GetPoints();
    int num_comp_points = component->GetNumberOfPoints();
    int num_comp_cells = component->GetNumberOfCells();

    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;
    for (int i = 0; i < num_comp_points; i++) {
      double pt[3];
      comp_points->GetPoint(i, pt);
      cx += pt[0];
      cy += pt[1];
      cz += pt[2];
    }

    center[0] = cx / num_comp_points;
    center[1] = cy / num_comp_points;
    center[2] = cz / num_comp_points;
    double d = (center[0]-position[0])*(center[0]-position[0]) +
         (center[1]-position[1])*(center[1]-position[1]) +
         (center[2]-position[2])*(center[2]-position[2]);

    if (d < min_d) {
      min_d = d;
      min_comp = component;
    }
    rid += 1;
  }
  /*
  while (true) {
    std::cout << "hi" << std::endl << std::flush;
    conn_filter->AddSpecifiedRegion(rid);
    conn_filter->Update();
    std::cout << "mmm" << std::endl << std::flush;
    auto component = vtkPolyData::New();
    std::cout << ".hi" << std::endl << std::flush;
    component->DeepCopy(conn_filter->GetOutput());
    std::cout << "..hi" << std::endl << std::flush;
    if (component->GetNumberOfCells() <= 0) {
      break;
    }
    conn_filter->DeleteSpecifiedRegion(rid);
    std::cout << "hi2" << std::endl << std::flush;
    auto clean_filter = vtkSmartPointer<vtkCleanPolyData>::New();
    clean_filter->SetInputData(component);
    clean_filter->Update();
    component = clean_filter->GetOutput();
    std::cout << "hi3" << std::endl << std::flush;

    auto comp_points = component->GetPoints();
    int num_comp_points = component->GetNumberOfPoints();
    int num_comp_cells = component->GetNumberOfCells();
    std::cout << "h4" << std::endl << std::flush;

    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;
    for (int i = 0; i < num_comp_points; i++) {
      double pt[3];
      comp_points->GetPoint(i, pt);
      cx += pt[0];
      cy += pt[1];
      cz += pt[2];
    }
    std::cout << "hi5" << std::endl << std::flush;

    center[0] = cx / num_comp_points;
    center[1] = cy / num_comp_points;
    center[2] = cz / num_comp_points;
    double d = (center[0]-position[0])*(center[0]-position[0]) +
         (center[1]-position[1])*(center[1]-position[1]) +
         (center[2]-position[2])*(center[2]-position[2]);

    if (d < min_d) {
      if (min_d > 1e5)
      {
  std::cout << "3" << std::endl << std::flush;
  // min_comp->Delete();
  std::cout << "3." << std::endl << std::flush;
      }
      min_d = d;
      min_comp = component;
    }
    else {
      std::cout << "4" << std::endl << std::flush;
      component->Delete();
      std::cout << "4." << std::endl << std::flush;
    }
    std::cout << "hi6" << std::endl << std::flush;
    rid += 1;
  }
  */
  return min_comp;
}

//-------------
// write_slice
//-------------
// Write a slice to a VTP file.
//
void Mesh::write_vtp(vtkPolyData* slice, std::string filename)
{
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(slice);
  writer->Update();
  writer->Write();
}

void Mesh::write_centerlines_and_fields(vtkPolyData* centerlines, std::string filename)
{
  for (auto it = flowrates_.begin(); it != flowrates_.end(); it++)
    centerlines->GetPointData()->AddArray(it->second);
  for (auto it = avg_pressures_.begin(); it != avg_pressures_.end(); it++)
    centerlines->GetPointData()->AddArray(it->second);
  centerlines->GetPointData()->AddArray(areas_);
  write_vtp(centerlines, filename);
}

