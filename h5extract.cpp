#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <H5Cpp.h>
#include <H5CommonFG.h>

using namespace std;
using namespace H5;

// dataspace from lengths and strides. Correct for the complex. strides must be >0
DataSpace dataspace_from_LS(int R, bool is_complex, hsize_t const *Ltot, hsize_t const *L, hsize_t const *S, hsize_t const *offset = NULL) {
   int rank = R + (is_complex ? 1 : 0);
   hsize_t totdimsf[rank], dimsf[rank], stridesf[rank], offsetf[rank]; 
   for (size_t u = 0; u < R; ++u) {
      offsetf[u] = (offset ? offset[u] : 0);
      dimsf[u] = L[u];
      totdimsf[u] = Ltot[u];
      stridesf[u] = S[u];
   }
   if (is_complex) {
      offsetf[rank - 1] = 0;
      dimsf[rank - 1] = 2;
      totdimsf[rank - 1] = 2;
      stridesf[rank - 1] = 1;
   }

   DataSpace ds = H5Screate_simple(rank, totdimsf, NULL);

   return ds;
}

// the dataspace corresponding to the array. Contiguous data only...
   DataSpace data_space_for_vector(std::vector<double> const &V) {
   hsize_t L[1], S[1];
   S[0] = 1;
   L[0] = V.size();
   return dataspace_from_LS(1, false, L, L, S);
}

   DataSpace data_space_for_vector(std::vector<complex<double>> const &V) {
   hsize_t L[1], S[1];
   S[0] = 1;
   L[0] = V.size();
   return dataspace_from_LS(1, true, L, L, S);
}

   void read_sigmak_into(vector<complex<double>>& V, string datasetPath, H5File& fp)
   {
    DataSet ds = fp.openDataSet(datasetPath.c_str());
    DataSpace d_space = ds.getSpace();

    static const unsigned int Rank = 3;
    hsize_t dims_out[Rank];
    d_space.getSimpleExtentDims(dims_out, NULL);
    V.resize(dims_out[0]*dims_out[1]);

    ds.read(&V[0], PredType::NATIVE_DOUBLE, data_space_for_vector(V), d_space, H5P_DEFAULT);
   }

   void read_vertex_into(vector<complex<double>>& V, string datasetPath, H5File& fp)
   {
    DataSet ds = fp.openDataSet(datasetPath.c_str());
    DataSpace d_space = ds.getSpace();

    static const unsigned int Rank = 4;
    hsize_t dims_out[Rank];
    d_space.getSimpleExtentDims(dims_out, NULL);
    V.resize(dims_out[0]*dims_out[1]*dims_out[2]);

    ds.read(&V[0], PredType::NATIVE_DOUBLE, data_space_for_vector(V), d_space, H5P_DEFAULT);
   }

   void read_into(vector<complex<double>>& V, string datasetPath, H5File& fp)
   {
    DataSet ds = fp.openDataSet(datasetPath.c_str());
    DataSpace d_space = ds.getSpace();

    static const unsigned int Rank = 2;
    hsize_t dims_out[Rank];
    d_space.getSimpleExtentDims(dims_out, NULL);
    V.resize(dims_out[0]);

    ds.read(&V[0], PredType::NATIVE_DOUBLE, data_space_for_vector(V), d_space, H5P_DEFAULT);
   }

   void read_into(vector<double>& V, string datasetPath, H5File& fp)
   {
    DataSet ds = fp.openDataSet(datasetPath.c_str());
    DataSpace d_space = ds.getSpace();

    static const unsigned int Rank = 1;
    hsize_t dims_out[Rank];
    d_space.getSimpleExtentDims(dims_out, NULL);
    V.resize(dims_out[0]);

    ds.read(&V[0], PredType::NATIVE_DOUBLE, data_space_for_vector(V), d_space, H5P_DEFAULT);
   }

int main(int argc, char* argv[])
{
    if (argc < 3) 
    {
      cout << "Usage: h5extract HDFFILE DATASET [BOSONIC FREQ FOR VERTEX]" << endl;
      exit(1);
    }

    string ifn = argv[1];
    H5File fp(ifn.c_str(),H5F_ACC_RDONLY);

    string datasetPath = argv[2];
   
    // work out the number of frequency indices in the grid
    Group dsgrids = fp.openGroup(datasetPath+"/grids");
    int num_freqs = dsgrids.getNumObjs();

   if (num_freqs == 3) // must be a vertex function
   {
    if (argc != 4) 
    {
      cout << "Need to specify bosonic frequency for vertex function" << endl;
      exit(1);
    }

    int bosonic_index = atoi(argv[3]);

    string datasetPath1 = datasetPath + "/data";
    string datasetPath2 = datasetPath + "/grids/0/values";
    string datasetPath3 = datasetPath + "/grids/1/values";
    string datasetPath4 = datasetPath + "/grids/2/values";
    vector<complex<double>> V;
    vector<complex<double>> grid0;
    vector<complex<double>> grid1;
    vector<complex<double>> grid2;
    read_vertex_into(V, datasetPath1, fp);
    read_into(grid0, datasetPath2, fp);
    read_into(grid1, datasetPath3, fp);
    read_into(grid2, datasetPath4, fp);

    int start = ((grid0.size()-1)/2 + bosonic_index)*grid1.size()*grid2.size();
    int end = ((grid0.size()-1)/2 + bosonic_index + 1)*grid1.size()*grid2.size();
    std::cout << start << " " << end << " " << grid0.size() << " " << grid1.size() << " " << grid2.size() << std::endl;

    for(int i = start; i < end; ++i)
    {
      int n1 = (i-start) % int(grid1.size());
      int n2 = (i-start) / int(grid1.size());
      std::cout << imag(grid1[n1]) << " " << imag(grid2[n2]) << " " << real(V[i]) << " " << imag(V[i]) << std::endl;

      if (n1 == grid1.size()-1)
         std::cout << std::endl;
    }
   }
   else if (num_freqs == 2) // must be a sigma(k,k')
   {
         string datasetPath1 = datasetPath + "/data";
         string datasetPath2 = datasetPath + "/grids/0/values";
         string datasetPath3 = datasetPath + "/grids/1/values";
         vector<complex<double>> V;
         vector<double> grid1;
         vector<double> grid2;
         read_sigmak_into(V, datasetPath1, fp);
         read_into(grid1, datasetPath2, fp);
         read_into(grid2, datasetPath3, fp);

         for(int i = 0; i < V.size(); ++i)
         {
            int n1 = i % int(grid1.size());
            int n2 = i / int(grid1.size());
            std::cout << grid1[n1] << " " << grid2[n2] << " " << real(V[i]) << " " << imag(V[i]) << std::endl;

            if (n1 == grid1.size()-1)
               std::cout << std::endl;
         }
      } 
      else if (num_freqs == 1) //must be a Green function
      {
         string datasetPath1 = datasetPath + "/data";
         string datasetPath2 = datasetPath + "/grids/0/values";
         vector<complex<double>> V;
         vector<complex<double>> grid;
         read_into(V, datasetPath1, fp);
         read_into(grid, datasetPath2, fp);

         for(int i = 0; i < V.size(); ++i)
            std::cout << imag(grid[i]) << " " << real(V[i]) << " " << imag(V[i]) << std::endl;
      }

   // close the HDF5 file
    fp.close();

return 0;
}
