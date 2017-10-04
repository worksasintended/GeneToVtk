#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "definitions.h"

using namespace std;
vector<Surface> readFile(string& filename, int& slice_points, int& surface_points, int& surfaceNumber);
void writeXML(string& outfile, int& timeSteps);
void print_memory_consumption( std::vector<Surface>& surfaces );
vector<float> find_limits(Surface& surface);
void interpolate_points(vector<SmallSurface>& volume, float& resolutionGoal);
void interpolate_slices(vector<SmallSurface>& s, float& resolutionGoal);
void interpolate_surfaces( std::vector<SmallSurface>& ss, float& resolutionGoal);
void map_to_grid(vector<SmallSurface>& surfaces, vector<float>& vtk_grid, vector<int>& ctr_grid, vector<float>& limits, int& x_gridResolution, int& y_gridResolution, int& z_gridResolution);
void writeVTK(string& outname, vector<float>& vtk_grid, vector<int>& ctr_grid, int& timeStep, int& x_gridResolution, int& y_gridResolution, int& z_gridResolution);


int main(int argc, char** argv){

//read argvs and check if argc is correct
  	if(argc != 8){
    		cout << "usage: ./bin_preprocessor <inputFile> <numberOfPointsInSurface> <numberOfPointsInSlice> <numberOfSurfaces> <gridResolution> <timeSteps> <outputFile> " << std::endl;
    		exit(-1);
  	}
	string filename = argv[1];
  	int surface_points;
  	int slice_points;
  	int surfaceNumber;
  	int grid_resolution;
	int timeSteps;

	try{
	  surface_points = stoi(argv[2]);
	  slice_points = stoi(argv[3]);
	  surfaceNumber = stoi(argv[4]);
	  grid_resolution = stoi(argv[5]);
	  timeSteps = stoi(argv[6]);
	}catch( ... ){
	  cout << "usage: ./bin_preprocessor <inputFile> <numberOfPointsInSurface> <numberOfPointsInSlice> <numberOfSurfaces> <gridResolution> <timeSteps> <outputFile> " << std::endl;
	  exit(-1);
	}
  	string outname = argv[7];



//read input file 
	vector<Surface> surfaces = readFile(filename, slice_points, surface_points, surfaceNumber);
	if(timeSteps > surfaces[0].slices[0].points[0].time_values.size()){
	    cout << "you asked for more time steps than the input file delivers" << endl;
	    cout << "usage: ./bin_preprocessor <inputFile> <numberOfPointsInSurface> <numberOfPointsInSlice> <numberOfSurfaces> <gridResolution> <timeSteps> <outputFile> " << std::endl;
	    exit(-1);
	}
#if 0
	ofstream out(outname);
	for(auto&& surface:surfaces){
	  for(auto&& slice:surface.slices){
	    for(auto&& point:slice.points){
	      out << point.x << "," << point.y << "," << point.z << "," << point.time_values[3] << endl;
	    }
	  }
	}
#endif
#if 1
//find limits in outer surface (sorted: x_min, x_max, y_min, ...., z_max)
	vector<float> limits = find_limits(surfaces[surfaceNumber-1]);
 //calculate resolution goal and grid dimensions
        float resolutionGoal = (limits[1]-limits[0])/grid_resolution/1.4;
	cout << "resolution goal: " << resolutionGoal << endl;
        int x_gridResolution = grid_resolution;
        int y_gridResolution = x_gridResolution*(limits[3]-limits[2])/(limits[1]-limits[0]);
        int z_gridResolution = x_gridResolution*(limits[5]-limits[4])/(limits[1]-limits[0]);
        cout << "dimensions of grid: " << x_gridResolution << "," << y_gridResolution << "," << z_gridResolution << endl;
	int memorySize = x_gridResolution*y_gridResolution*z_gridResolution*(8)/1024/1024;
	cout << "memory consumption of grid: " << memorySize << " MiB" << endl;

// call interpolation routine for each timestep
	for(int timeStep=0; timeStep < timeSteps; ++timeStep){
	  cout << "time step: " << timeStep << endl;
	  //create output grid and counter grid
	  vector<float> vtk_grid((x_gridResolution+1)*(y_gridResolution+1)*(z_gridResolution+1));
	  fill(vtk_grid.begin(), vtk_grid.end(), 0.);
	  vector<int> ctr_grid((x_gridResolution+1)*(y_gridResolution+1)*(z_gridResolution+1));
	  fill(ctr_grid.begin(), ctr_grid.end(), 0);
	  //take a small volume part between two slices
	  for (int startSlice=0; startSlice < surfaces[0].slices.size()-1; ++startSlice ){

	    vector<SmallSurface> newVolume;
	    for(auto& surface : surfaces){
	      SmallSurface newSurface;
	      for(int i=0; i<2 ; ++i){
		Slice tempSlice= surface.slices[startSlice+i];
		SmallSlice tempSmallSlice;
		for(int j=0; j  <tempSlice.points.size(); ++j){
		  SmallPoint tempSmallPoint;
		  tempSmallPoint.x=tempSlice.points[j].x;
		  tempSmallPoint.y=tempSlice.points[j].y;
		  tempSmallPoint.z=tempSlice.points[j].z;
		  tempSmallPoint.value=tempSlice.points[j].time_values[timeStep];
		  tempSmallSlice.points.push_back(tempSmallPoint);
		}
		newSurface.slices.push_back(tempSmallSlice);
	      }
	      newVolume.push_back(newSurface);
	    }
	    //remove all timesteps despite the one used
	    interpolate_points(newVolume, resolutionGoal);
	    interpolate_slices(newVolume, resolutionGoal);
	    interpolate_surfaces(newVolume, resolutionGoal);
	    cout << "done volume part interpolation " << startSlice << " / " << surfaces[0].slices.size()-1 << endl;
	//    cout << "test" << endl;
	    //print_memory_consumption(newVolume);
	    map_to_grid(newVolume, vtk_grid, ctr_grid, limits, x_gridResolution, y_gridResolution, z_gridResolution);  
	  }

	  
	  writeVTK(outname, vtk_grid, ctr_grid, timeStep, x_gridResolution, y_gridResolution, z_gridResolution);  
	}

//write xml file
	writeXML(outname, timeSteps);
#endif


	return 0;
}
void writeVTK(string& outname, vector<float>& vtk_grid, vector<int>& ctr_grid, int& timeStep, int& x_gridResolution, int& y_gridResolution, int& z_gridResolution){
  stringstream step;
  step << setfill('0') << setw(4) << timeStep;
  string outputVtk = outname + step.str() + ".vtk";
  ofstream out(outputVtk);
  out << "# vtk DataFile Version 2.0" << endl << "interpolated data" << endl << "ASCII"<< endl << "DATASET STRUCTURED_POINTS" << endl << "DIMENSIONS " << x_gridResolution << " " << y_gridResolution << " " <<  z_gridResolution << endl << "SPACING 1 1 1" << endl << "ORIGIN 0 0 0 "<< endl << "POINT_DATA " << x_gridResolution*y_gridResolution*z_gridResolution << endl << "SCALARS simulatedValues float" << endl << "LOOKUP_TABLE default" << endl;

   for(int i=0; i<z_gridResolution; ++i){
    for(int j=0; j< y_gridResolution; ++j){
      for(int k=0; k< x_gridResolution; ++k){
	 double temp_value;
	 if(ctr_grid[k+(x_gridResolution)*j+(x_gridResolution)*(y_gridResolution)*i]==0){ temp_value=0;}
	 else{temp_value=vtk_grid[k+(x_gridResolution)*j+(x_gridResolution)*(y_gridResolution)*i]/ctr_grid[k+(x_gridResolution)*j+(x_gridResolution)*(y_gridResolution)*i];}
	 out << temp_value << " "; 
      }
    }
    out << endl;
  }
  cout << "file " << outputVtk << " written" << endl;
 

}


void map_to_grid(vector<SmallSurface>& surfaces, vector<float>& vtk_grid, vector<int>& ctr_grid, vector<float>& limits, int& x_gridResolution, int& y_gridResolution, int& z_gridResolution){
  for( auto&& surface : surfaces){
    for( auto&& slice : surface.slices){
      for(auto&& point: slice.points){
	int x_grid_position=(point.x-limits[0])/(limits[1]-limits[0])*x_gridResolution;
	int y_grid_position=(point.y-limits[2])/(limits[3]-limits[2])*y_gridResolution;
	int z_grid_position=(point.z-limits[4])/(limits[5]-limits[4])*z_gridResolution;
	//cout << x_grid_position << " " << y_grid_position << " " << z_grid_position << " " << point.z << endl;
	float value = point.value;
//	cout << point.time_values.size()<< " " ;
        vtk_grid[x_grid_position+x_gridResolution*y_grid_position+x_gridResolution*y_gridResolution*z_grid_position]+=value;
	ctr_grid[x_grid_position+x_gridResolution*y_grid_position+x_gridResolution*y_gridResolution*z_grid_position]++;

	}
      }
  } 
  
}

void interpolate_surfaces( std::vector<SmallSurface>& ss, float& resolutionGoal){
    float resolution = sqrt(pow(ss[0].slices[0].points[0].x-ss[1].slices[0].points[0].x, 2)+pow(ss[0].slices[0].points[0].y-ss[1].slices[0].points[0].y, 2) + pow(ss[0].slices[0].points[0].z-ss[1].slices[0].points[0].z, 2));
    //cout << "first resolution guess: " << resolution << endl;
    for(int i=0; i<ss.size()-1;++i){
        SmallSurface& s1=ss[i];
        SmallSurface& s2=ss[i+1];
        for(int j=0; j<ss[0].slices.size(); ++j){
            SmallSlice& sl1=s1.slices[j];
            SmallSlice& sl2=s2.slices[j];
            for(int k=0; k<ss[0].slices[0].points.size(); ++k){
                float newResolution=sqrt(pow(sl1.points[k].x-sl2.points[k].x,2)+pow(sl1.points[k].y-sl2.points[k].y,2)+pow(sl1.points[k].z-sl2.points[k].z,2));
                if(newResolution > resolution){
                    resolution=newResolution;
                }

            }
        }
    }
  //  cout << "final resolution: " << resolution << endl;

    while(resolution > resolutionGoal){
    std::vector<SmallSurface> new_volume = ss;
    for (int i = 0; i <ss.size()-1; ++i){
      SmallSurface new_surface;
      auto& surface_a=ss[i];
      auto& surface_b=ss[i+1];
     // std::cout << i << std::endl;
      for (int j = 0; j < surface_a.slices.size(); j++){
	SmallSlice new_slice;
	auto& slice_a = surface_a.slices[j];
	auto& slice_b = surface_b.slices[j];
	//std::cout << j << endl;
	for (int k=0; k< slice_a.points.size(); k ++){
	  auto& point_a = slice_a.points[k];
	  auto& point_b = slice_b.points[k];
	  SmallPoint new_point;
	  new_point.x = (point_a.x + point_b.x) / 2;
	  new_point.y = (point_a.y + point_b.y) / 2;
	  new_point.z = (point_a.z + point_b.z) / 2;
	  new_point.value= (point_a.value + point_b.value) / 2;	  
	  new_slice.points.push_back( new_point);
	}

	new_surface.slices.push_back( new_slice);
      }
	    
      new_volume.insert( new_volume.begin()+ 2.*i+1,new_surface);
    }

    ss= new_volume;
    resolution = resolution /2 ;
    }
}

void interpolate_slices( vector<SmallSurface>& surfaces, float& resolutionGoal)
{
  float resolution = sqrt(pow(surfaces[0].slices[0].points[0].x-surfaces[0].slices[1].points[0].x, 2)+pow(surfaces[0].slices[0].points[0].y-surfaces[0].slices[1].points[0].y, 2) + pow(surfaces[0].slices[0].points[0].z-surfaces[0].slices[1].points[0].z, 2));
  for(auto&& surface:surfaces){
    for(int i=0; i<surface.slices.size()-1; ++i){
      SmallSlice& s1=surface.slices[i];
      SmallSlice& s2=surface.slices[i+1];
      for(int j=0; j<surface.slices[0].points.size(); ++j){
	SmallPoint& p1=s1.points[j];
	SmallPoint& p2=s2.points[j];
	float newResolution = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
	if(newResolution > resolution){
	  resolution=newResolution;
	}
      }
    }

  }
#if 1
  while(resolution > resolutionGoal){
      for(auto& s : surfaces){
	//std::cout << "slices = " << s.slices.size() << std::endl;
	SmallSurface new_surface = s;
  	for (int i = 0; i < s.slices.size()-1; ++i){
    		auto& slice_a = s.slices[i];
    		auto& slice_b = s.slices[i+1];
   		SmallSlice new_slice;
    		for (int pid = 0; pid < slice_a.points.size(); ++pid){
      			auto& point_a = slice_a.points[pid];
      			auto& point_b = slice_b.points[pid];
			SmallPoint new_point;
      			new_point.x = (point_a.x + point_b.x) / 2;
      			new_point.y = (point_a.y + point_b.y) / 2;
      			new_point.z = (point_a.z + point_b.z) / 2;
			new_point.value= (point_a.value + point_b.value) / 2;	  
      			new_slice.points.push_back( new_point );
    		}
    		new_surface.slices.insert( new_surface.slices.begin() + i*2 + 1, new_slice );
  	}
    
	//  std::cout << "slices = " << new_surface.slices.size() << std::endl;
  	s = new_surface;
      }
    resolution = resolution /2 ;
  }
#endif

}

void interpolate_points(vector<SmallSurface>& surfaces, float& resolutionGoal){
  vector<SmallSurface> newVolume = surfaces;
  float resolution = sqrt(pow(surfaces[0].slices[0].points[0].x-surfaces[0].slices[0].points[1].x, 2)+pow(surfaces[0].slices[0].points[0].y-surfaces[0].slices[0].points[1].y, 2) + pow(surfaces[0].slices[0].points[0].z-surfaces[0].slices[0].points[1].z, 2));
  for(auto&& surface:surfaces){
    for(auto&& slice:surface.slices){
      for(int i=0; i<slice.points.size()-1; ++i){
	SmallPoint& p1=slice.points[i];
	SmallPoint& p2=slice.points[i+1];
	float newResolution = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
	if(newResolution > resolution){
	  resolution=newResolution;
	}
      }
    }
  }


  while(resolution > resolutionGoal){
    for (int i = 0; i < surfaces.size(); ++i){
      for(int j=0; j< surfaces[0].slices.size(); ++j){
	for(int k=0; k< surfaces[0].slices[0].points.size()-1; ++k){
	  SmallPoint& pointA=surfaces[i].slices[j].points[k];
	  SmallPoint& pointB=surfaces[i].slices[j].points[k+1];
	  SmallPoint newPoint;
	  newPoint.x = (pointA.x + pointB.x)/2.;
	  newPoint.y = (pointA.y + pointB.y)/2.;
	  newPoint.z = (pointA.z + pointB.z)/2.;
	  newPoint.value = (pointA.value + pointB.value) / 2;
	  newVolume[i].slices[j].points.insert(newVolume[i].slices[j].points.begin() + k*2 +1, newPoint);

	}
      }
    }
    surfaces = newVolume;
    resolution = resolution /2;
  }

}

vector<Surface> readFile(string& filename, int& slice_points, int& surface_points, int& surfaceNumber){
	ifstream in(filename);
	if(!in.good()){
		cout << "input file cannot be read" << endl << "usage: ./bin_preprocessor <inputFile> <numberOfPointsInSurface> <numberOfPointsInSlice> <gridResolution> <timeSteps> <outputFile> " << std::endl; 
		exit(-1);
	}
	vector<Surface> surfaces;
	//skip the first lines that are comments
	for(int i=0; i<5; ++i){
		string trash;
		getline(in, trash);
		cout << trash << endl;
	}
	for(int i=0; i< surfaceNumber ; ++i ) {
		std::string line;
    		Surface surface;
    		for (int i = 0; i < surface_points/slice_points; ++i){
      			Slice s;
      			for (int j = 0; j < slice_points; ++j){
				if( getline(in, line) ){
	  				stringstream lstr( line ) ;
	  				Point p;
	 				lstr >> p.phi >> p.theta >> p.x >> p.y >> p.z;
	  				float value;
					while( lstr >> value ) {
	    					p.time_values.push_back(value);
	 				}
	  				s.points.push_back( p );
				} 
      			}
      			surface.slices.push_back( s );
    		}
		surfaces.push_back( surface );
  	}
  	//seems like there is an additional empty entry at the end --> remove it
  	//surfaces.pop_back();
	//surfaces.pop_back();
	cout << "input file read" << endl;
	print_memory_consumption(surfaces); 
	return surfaces; 
}

void writeXML(string& outfile, int& timeSteps){
// get number of timesteps 
	string filename=outfile+".xml";
	ofstream out(filename);
	string tempString;
	out << "<?xml version=\"1.0\"?>" << endl << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" <<  endl << " <Collection>" << endl;
	for(int i=0; i<timeSteps; ++i){
		stringstream step; 
		step << setfill('0') << setw(4) << i;
		string tempString = outfile + step.str() +".vtk";
		out << "  <DataSet timestep=\"" << i <<"\" group=\"\" part=\"0\" file=\""<< tempString << "\"/> " << endl;
	}
	out << " </Collection>" << endl << "</VTKFile>";
	cout << "output files: " << filename << " and corresponding .vtk\'s for "<< timeSteps << " time steps" << endl; 
}

void print_memory_consumption( std::vector<Surface>& surfaces ){
	long ctr = 0;

	for( auto& surface : surfaces ) {
      		for( auto& slice : surface.slices ) {
	  		for( auto& point : slice.points ) {
	    			ctr += (5 + point.time_values.size()) * 4;
	  		}
      		}
    	}
    	std::cout << "memory consumption " << ctr/1024/1024 << " MiB" << endl;
}

vector<float> find_limits(Surface& surface){
	vector<float> limits;;
	float limit = surface.slices[0].points[0].x;
	for(auto&& slice:surface.slices){
	  for(auto&& point:slice.points){
	    if(point.x < limit){limit=point.x;}
	  } 
	}
	limits.push_back(limit);
	limit = surface.slices[0].points[0].x;
	for(auto&& slice:surface.slices){
	  for(auto&& point:slice.points){
	    if(point.x > limit){limit=point.x;}
	  } 
	}
	limits.push_back(limit);
	limit = surface.slices[0].points[0].y;
	for(auto&& slice:surface.slices){
	  for(auto&& point:slice.points){
	    if(point.y < limit){limit=point.y;}
	  } 
	}
	limits.push_back(limit);
	limit = surface.slices[0].points[0].y;
	for(auto&& slice:surface.slices){
	  for(auto&& point:slice.points){
	    if(point.y > limit){limit=point.y;}
	  } 
	}
	limits.push_back(limit);
	limit = surface.slices[0].points[0].z;
	for(auto&& slice:surface.slices){
	  for(auto&& point:slice.points){
	    if(point.z < limit){limit=point.z;}
	  } 
	}
	limits.push_back(limit);
	limit = surface.slices[0].points[0].z;
	for(auto&& slice:surface.slices){
	  for(auto&& point:slice.points){
	    if(point.z > limit){limit=point.z;}
	  } 
	}
	limits.push_back(limit);
	cout << "limits are: " << "	x:[" << limits[0] << "," << limits[1] << "], y:[" << limits[2] << "," << limits[3] << "], z:[" << limits[4] << "," << limits[5] << "]" << endl;
	return limits;

}
