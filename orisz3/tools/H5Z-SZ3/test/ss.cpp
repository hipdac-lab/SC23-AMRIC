#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <dlfcn.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "hdf5.h"
#include "H5Z_SZ3.hpp"
#include <cstring>
#define DATASET "/level_0/data:datatype=0"
#define MAX_CHUNK_SIZE 4294967295 //2^32-1
#define BSIZE 16
#define LEVEL 0
#define XSIZE 128
#define YSIZE 128
#define ZSIZE 128

bool boxesInsec(int* box1, int* box2) {
        if (box1[3] < box2[0] || box2[3] < box1[0]) {
        return false;
    }

    // Check y dimension
    if (box1[4] < box2[1] || box2[4] < box1[1]) {
        return false;
    }

    // Check z dimension
    if (box1[5] < box2[2] || box2[5] < box1[2]) {
        return false;
    }

    // Boxes overlap in all three dimensions
    return true;
}

int main(int argc, char * argv[])
{
        //dcdc test
        // const char *chunk_env = getenv("HDF5_CHUNK_SIZE");
        // hsize_t chunk = atoi(chunk_env);

        // printf("chunk size = %zu\n", chunk);

        int dimSize = 0;
        size_t r5=0,r4=0,r3=0,r2=0,r1=0;
        long nbEle = 0;
        char hdf5FilePath[640], outputFilePath[640];
        hid_t file, dset, dcpl, space_id, dtype; /*Handles*/
        H5Z_filter_t filter_id = 0;
        herr_t status;
        H5T_class_t type_class;
        H5T_sign_t dsign;
        H5T_order_t dorder;

        htri_t avail;
        char filter_name[80];
        unsigned int flags = 0;
        size_t nelmts = 0, dsize;
        unsigned int values_out[7] = {0,0,0,0,0,0,0}; //at most 7 parameters

        if(argc < 2)
        {
                printf("Test case: dszFromHDF5 [hdf5FilePath]\n");
                printf("Example 1: dszFromHDF5 testdata/x86/testfloat_8_8_128.dat.sz.hdf5\n");
                printf("Example 2: dszFromHDF5 testdata/x86/testint32_8x8x8.dat.sz.hdf5\n");
                exit(0);
        }

        sprintf(hdf5FilePath, "%s", argv[1]);
        sprintf(outputFilePath, "%s.out.h5", hdf5FilePath);

        /*Open the hdf5 file with SZ-compressed data*/
    file = H5Fopen(hdf5FilePath, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen(file, DATASET, H5P_DEFAULT);

    /*Retrieve boxDataset creation property list.*/
    dcpl = H5Dget_create_plist(dset);

    /*Check that filter is not registered with the library yet*/


        space_id = H5Dget_space(dset);
        nbEle = H5Sget_simple_extent_npoints(space_id);
        std::cout << "nbEle: " << nbEle << std::endl;

        /*partio I/O*/
        // hsize_t count[1], offsets[1];
        // count[0] = nbEle/NFIELD;
        // offsets[0] = MYFIELD * nbEle/NFIELD;
        // std::cout << "offsets[0] " << offsets[0] << std::endl;
        // status = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offsets, NULL, count, NULL);
        // hid_t memspace_id = H5Screate_simple(1, count, NULL);



        if((dtype = H5Dget_type(dset)) < 0)
                printf("Error: H5Dget_type(dset) < 0\n");

        /*Read the data using the default properties.*/
        printf("....Reading SZ compressed data .....................\n");

        if((type_class = H5Tget_class(dtype)) < 0)
        {
                printf("Error: H5Tget_class<0\n");
                exit(0);
        }
        if (0 == (dsize = H5Tget_size(dtype)))
        {
                printf("Error: H5Tget_size==0\n");
                exit(0);
        }

        if((dorder = H5Tget_order(dtype)) < 0)
                printf("Error: H5Tget_order<0\n");




        double* realData = (double*)malloc(sizeof(double)*nbEle);

        switch (type_class)
        {
        case H5T_FLOAT:
                if (H5Tequal(dtype, H5T_IEEE_F32BE) == 1 || H5Tequal(dtype, H5T_IEEE_F32LE) == 1
                || H5Tequal(dtype, H5T_NATIVE_FLOAT) == 1)
                {
                        printf("data type: float\n");
                        float* data = (float*)malloc(sizeof(float)*nbEle);
                        if(dorder==H5T_ORDER_LE)
                                status = H5Dread(dset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                        else //H5T_ORDER_BE
                                status = H5Dread(dset, H5T_IEEE_F32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                        /*Print the first 20 data values to check the correctness.*/
                        int i;
                        printf("reconstructed data = ");
                        for(i=0;i<20;i++)
                                printf("%f ", data[i]);
                        printf("\n");
                        free(data);
        }
                else //64bit: double
                {
                        printf("data type: double\n");
                        double* data = (double*)malloc(sizeof(double)*nbEle);
                        if(dorder==H5T_ORDER_LE)
                                status = H5Dread(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                // status = H5Dread(dset, H5T_IEEE_F64LE, memspace_id, space_id, H5P_DEFAULT, data);
                        else
                                status = H5Dread(dset, H5T_IEEE_F64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                // status = H5Dread(dset, H5T_IEEE_F64LE, memspace_id, space_id, H5P_DEFAULT, data);

                        for(auto i = 0; i < nbEle; ++i) {
                                realData[i]=data[i];
                        }




                        // std::ofstream fout("test.bin", std::ios::binary);
                        // fout.write(reinterpret_cast<const char*>(realData), realSum * sizeof(double));
                        // fout.close();

                        // for(int i=0;i<realSum;i++){
                        //      if (realData[i]<1100)
                        //              std::cout << i << std::endl;
                        // }

                        /*Print the first 10 data values to check the correctness.*/
                        // int i;
                        // printf("reconstructed data = ");
                        // for(i=0;i<20;i++)
                        //      printf("%f ", data[i]);
                        // printf("\n");
                        // free(data);
                        double sum = 0;
                        for(int i=0;i<nbEle/10;i++)
                                sum+=data[i];
                        std::cout << sum << std::endl;
                }
                break;
        case H5T_INTEGER:
                if (0 > (dsign = H5Tget_sign(dtype)))
                {
                        printf("Error in calling H5Tget_sign(type_id)....\n");
                        exit(0);
                }
                if(dsign == H5T_SGN_NONE) //unsigned
                {
                        if(dsize==1)
                        {
                                printf("data type: unsigned char\n");
                                unsigned char* data = (unsigned char*)malloc(sizeof(unsigned char)*nbEle);
                                if(dorder==H5T_ORDER_LE)
                                        status = H5Dread(dset, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                else
                                        status = H5Dread(dset, H5T_STD_U8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                int i;
                                printf("reconstructed data = ");
                                for(i=0;i<20;i++)
                                        printf("%d ", data[i]);
                                printf("\n");
                                free(data);
                        }
                        else if(dsize==2)
                        {
                                printf("data type: unsigned short\n");
                                unsigned short* data = (unsigned short*)malloc(sizeof(unsigned short)*nbEle);
                                if(dorder==H5T_ORDER_LE)
                                        status = H5Dread(dset, H5T_STD_U16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                else
                                        status = H5Dread(dset, H5T_STD_U16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                int i;
                                printf("reconstructed data = ");
                                for(i=0;i<20;i++)
                                        printf("%d ", data[i]);
                                printf("\n");
                                free(data);
                        }
                        else if(dsize==4)
                        {
                                printf("data type: unsigned int\n");
                                unsigned int* data = (unsigned int*)malloc(sizeof(unsigned int)*nbEle);
                                if(dorder==H5T_ORDER_LE)
                                        status = H5Dread(dset, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                else
                                        status = H5Dread(dset, H5T_STD_U32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                int i;
                                printf("reconstructed data = ");
                                for(i=0;i<20;i++)
                                        printf("%d ", data[i]);
                                printf("\n");
                                free(data);
                        }
                        else if(dsize==8)
                        {
                                printf("data type: unsigned long\n");
                                unsigned long* data = (unsigned long*)malloc(sizeof(unsigned long)*nbEle);
                                if(dorder==H5T_ORDER_LE)
                                        status = H5Dread(dset, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                else
                                        status = H5Dread(dset, H5T_STD_U64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                int i;
                                printf("reconstructed data = ");
                                for(i=0;i<20;i++)
                                        printf("%ld ", data[i]);
                                printf("\n");
                                free(data);
                        }
                }
                else
                {
                        if(dsize==1)
                        {
                                printf("data type: char\n");
                                char *data = (char*)malloc(sizeof(char)*nbEle);
                                if(dorder==H5T_ORDER_LE)
                                        status = H5Dread(dset, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                else
                                        status = H5Dread(dset, H5T_STD_I8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                int i;
                                printf("reconstructed data = ");
                                for(i=0;i<20;i++)
                                        printf("%d ", data[i]);
                                printf("\n");
                                free(data);
                        }
                        else if(dsize==2)
                        {
                                printf("data type: short\n");
                                short *data = (short*)malloc(sizeof(short)*nbEle);
                                if(dorder==H5T_ORDER_LE)
                                        status = H5Dread(dset, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                else
                                        status = H5Dread(dset, H5T_STD_I16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                int i;
                                printf("reconstructed data = ");
                                for(i=0;i<20;i++)
                                        printf("%d ", data[i]);
                                printf("\n");
                                free(data);
                        }
                        else if(dsize==4)
                        {
                                printf("data type: int\n");
                                int *data = (int*)malloc(sizeof(int)*nbEle);
                                if(dorder==H5T_ORDER_LE)
                                        status = H5Dread(dset, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                else
                                        status = H5Dread(dset, H5T_STD_I32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                int i;
                                printf("reconstructed data = ");
                                for(i=0;i<20;i++)
                                        printf("%d ", data[i]);
                                printf("\n");
                                free(data);
                        }
                        else if(dsize==8)
                        {
                                printf("data type: long\n");
                                long *data = (long*)malloc(sizeof(long)*nbEle);
                                if(dorder==H5T_ORDER_LE)
                                        status = H5Dread(dset, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                else
                                        status = H5Dread(dset, H5T_STD_I64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                                int i;
                                printf("reconstructed data = ");
                                for(i=0;i<20;i++)
                                        printf("%ld ", data[i]);
                                printf("\n");
                                free(data);
                        }
                }

                break;
        default:
                printf("Error: H5Z-SZ supports only float, double or integers.\n");
                exit(0);
        }


        hid_t boxDataset_1, boxDatatype_1, boxDataspace_1;
    herr_t boxStatus_1;

    boxDataset_1 = H5Dopen2(file, "/level_1/boxes", H5P_DEFAULT);
    boxDatatype_1 = H5Dget_type(boxDataset_1);
        boxDataspace_1 = H5Dget_space(boxDataset_1);

        int boxLen_1 = H5Sget_simple_extent_npoints(boxDataspace_1);

        // std::cout << "boxLen_1: " << boxLen_1 << std::endl;
        int box_1[boxLen_1][6];

    boxStatus_1 = H5Dread(boxDataset_1, boxDatatype_1, H5S_ALL, H5S_ALL, H5P_DEFAULT, box_1);

        for (int j = 0; j < boxLen_1; ++j) {
                for (int i = 0; i < 6; ++i){
                        if (i < 3)
                                box_1[j][i] = box_1[j][i]/2;
                        else
                                box_1[j][i] = (box_1[j][i]+1)/2 - 1;
                        // std::cout << box_1[j][i] << std::endl;
                }
                // std::cout << " " << std::endl;
        }

    boxStatus_1 = H5Dclose(boxDataset_1);

        hid_t boxDataset_0, boxDatatype_0, boxDataspace_0;
    herr_t boxStatus_0;

    boxDataset_0 = H5Dopen2(file, "/level_0/boxes", H5P_DEFAULT);
    boxDatatype_0 = H5Dget_type(boxDataset_0);
        boxDataspace_0 = H5Dget_space(boxDataset_0);

        int boxLen_0 = H5Sget_simple_extent_npoints(boxDataspace_0);

        // std::cout << "boxLen_0: " << boxLen_0 << std::endl;
        int box_0[boxLen_0][6];

    boxStatus_0 = H5Dread(boxDataset_0, boxDatatype_0, H5S_ALL, H5S_ALL, H5P_DEFAULT, box_0);
        // for (int j = 0; j < boxLen_0; ++j) {
        //      for (int i = 0; i < 6; ++i){
        //              std::cout << box_0[j][i] << std::endl;
        //      }
        //      std::cout << " " << std::endl;
        // }
    boxStatus_0 = H5Dclose(boxDataset_0);

        /*dcdc get offset*/
        hid_t offDataset, offDatatype, offDataspace;
    herr_t offStatus;
    offDataset = H5Dopen2(file, "/level_0/data:offsets=0", H5P_DEFAULT);
    offDatatype = H5Dget_type(offDataset);
        offDataspace = H5Dget_space(offDataset);
        int offLen = H5Sget_simple_extent_npoints(offDataspace);
        // std::cout << "offLen: " << offLen << std::endl;
        long offset[offLen];
    offStatus = H5Dread(offDataset, offDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, offset);
        std::cout << "offset: ";
        for (int i = 0; i < offLen; ++i)
                std::cout <<  offset[i] << " ";
        std::cout << std::endl;
        offStatus = H5Dclose(offDataset);

        status = H5Pclose(dcpl);
        status = H5Dclose(dset);
        status = H5Fclose(file);

        double *testArr = new double[XSIZE*YSIZE*ZSIZE]{1.0};
        for (int i = 0; i < XSIZE*YSIZE*ZSIZE; ++i)
                testArr[i] = 1;

        for (int i = 0; i < boxLen_1; ++i) {
                        for (size_t z=box_1[i][2]; z<=box_1[i][5]; ++z)
                                for (size_t y=box_1[i][1]; y<=box_1[i][4]; ++y)
                                        for (size_t x=box_1[i][0]; x<=box_1[i][3]; ++x) {
                                                testArr[x + y*XSIZE + z*XSIZE*YSIZE] = 0;
                                        }
        }


        long cnt = 0;


        /*bs decomp*/
        for (int j = 0; j < boxLen_0; ++j) {
                for (size_t z=box_0[j][2]; z<=box_0[j][5]; ++z)
                        for (size_t y=box_0[j][1]; y<=box_0[j][4]; ++y)
                                for (size_t x=box_0[j][0]; x<=box_0[j][3]; ++x) {
                                        if(testArr[x + y*XSIZE + z*XSIZE*YSIZE] != 0) {
                                                //testArr[x + y*XSIZE + z*XSIZE*YSIZE]= realData[cnt];
                                                testArr[x + y*XSIZE + z*XSIZE*YSIZE] = realData[offset[j]+cnt];
                                        }
                                        cnt++;
                                }
                cnt = 0;
        }
        double nastSum=0;
        for (int i = 0; i < XSIZE*YSIZE*ZSIZE; ++i)
                nastSum+=testArr[i];
        std::cout << "Sum: " << nastSum << std::endl;


        std::cout << "cnt: " << cnt << std::endl;
        std::ofstream tout("ss.raw", std::ios::binary);
        tout.write(reinterpret_cast<const char*>(testArr), XSIZE*YSIZE*ZSIZE * sizeof(double));
        tout.close();

        free(realData);
        return 0;
}

