#include "output.h"
//determines the type based off of a user string
int determineType(char* type)
{
  if(strcmp(type,"double*") == 0)
  {
    return 0;
  }else if(strcmp(type,"double**") == 0)
  {
    return 1;
  }else if(strcmp(type,"double***") == 0)
  {
    return 2;
  }else if(strcmp(type,"double****") == 0)
  {
    return 3;
  }else if(strcmp(type,"int*") == 0)
  {
    return 4;
  }else if(strcmp(type,"int**") == 0)
  {
    return 5;
  }else if(strcmp(type,"int***") == 0)
  {
    return 6;
  }else if(strcmp(type,"int****") == 0)
  {
    return 7;
  }else
  {
    return -1;
  }
};
//function for writing a dataset to HDF5: arr is the file to be written, type denotes whether it is a double (0) or int (1) (can be expanded), store is the file being written to,
//name is the name of the dataset to be created for reference, dimnum denotes the number of dimensions in the set, dimlist is the size of those dimensions
void writeArr(void* arr, int type, char* filename, char* name, int dimnum, int* dimlist)
{
  hid_t file_id, set_id, dataspace_id;
  hsize_t* dims = (hsize_t*)dimlist;
  herr_t status;
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple(dimnum, dims, NULL);
  set_id = H5Dcreate2(file_id, name, H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(set_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
  status = H5Dclose(set_id);
  status = H5Sclose(dataspace_id);
  status = H5Fclose(file_id);
}