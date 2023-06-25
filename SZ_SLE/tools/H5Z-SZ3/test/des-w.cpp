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
#define ZSIZE 1024

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
        char meta_name[320];
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
        avail = H5Zfilter_avail(H5Z_FILTER_SZ3);
                printf("sz filter is available.\n");

        space_id = H5Dget_space(dset);
        nbEle = H5Sget_simple_extent_npoints(space_id);
        std::cout << "nbEle: " << nbEle << std::endl;

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

        int finest_level;
        FILE* fFile = fopen("../meta/f.txt", "r+");
        if (fFile == NULL) {
                printf("Could not open f file!\n");
                return 1;
        }
        fscanf(fFile, "%d", &finest_level);
        fclose(fFile);

        int ncomp;
        FILE* nFile = fopen("../meta/ncomp.txt", "r+");
        if (nFile == NULL) {
                printf("Could not open ncomp file!\n");
                return 1;
        }
        fscanf(nFile, "%d", &ncomp);
        fclose(nFile);

        int nProcs;
        FILE* pFile = fopen("../meta/p.txt", "r+");
        if (pFile == NULL) {
                printf("Could not open p file!\n");
                return 1;
        }
        fscanf(pFile, "%d", &nProcs);
        fclose(pFile);

        int realNprocs = 0;
        sprintf(meta_name, "../meta/realp_%d.txt", LEVEL);
        FILE* rpFile = fopen(meta_name, "r+");
        if (rpFile == NULL) {
                printf("Could not open p file!\n");
                return 1;
        }
        fscanf(rpFile, "%d", &realNprocs);
        fclose(rpFile);

        std::vector<long> procLength(nProcs);
        std::vector<long> procOffsets(nProcs);
        long realLength;



        if (LEVEL == finest_level)
        {
                sprintf(meta_name, "../meta/meta_%d.txt", LEVEL);
                FILE* fp = fopen(meta_name, "r+");
                if (fp == NULL) {
                        printf("Could not open file!\n");
                        return 1;
                }
                int i = 0;
                while (fscanf(fp, "%ld", &realLength) == 1){
                        procLength[i]=realLength;
                        ++i;
                }
                fclose(fp);
        } else {
                for (int i=0; i<nProcs; ++i) {
                        sprintf(meta_name, "../meta/meta_%d_%d.txt", LEVEL, i);
                        FILE* file = fopen(meta_name, "r+");
                        if (file == NULL) {
                                printf("Could not open file!\n");
                                return 1;
                        }
                        fscanf(file, "%ld", &realLength);
                        procLength[i]=realLength;
                        fclose(file);
                }
        }

        // procLength.erase(std::remove(procLength.begin(), procLength.end(), 0), procLength.end());
        std::cout << "procLength.length: " << procLength.size() << std::endl;
        std::cout << "realprocLength.length: " << realNprocs << std::endl;

        auto maxBuf = nbEle/(realNprocs*ncomp);
        std::cout << "maxBuf: " << maxBuf << std::endl;

        long long realNbEle(0);
        for(auto it = procLength.begin(); it != procLength.end(); ++it) {
                std::cout << *it << " ";
                realNbEle += *it;
        }
        realNbEle = realNbEle*ncomp;
        std::cout << std::endl;
        std::cout << "realNbEle: " << realNbEle << std::endl;
        double* realData = (double*)malloc(sizeof(double)*realNbEle);

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
                        else
                                status = H5Dread(dset, H5T_IEEE_F64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

                        long long writeDataSize = 0;
                        int cpCnt = 0;
                        for(auto i = 0; i < procLength.size(); ++i) {
                                if (procLength[i]>-1) {
                                        procOffsets[i] = writeDataSize;
                                        for(auto j = 0; j < ncomp; ++j) {
                                                // std::cout << "writeDataSize: " << (j+cpCnt*ncomp)*maxBuf << " " << writeDataSize << std::endl;
                                                memcpy(realData+writeDataSize, data+((j+cpCnt*ncomp)*maxBuf), procLength[i] * sizeof(double));
                                                writeDataSize += procLength[i];
                                        }
                                        cpCnt++;
                                }
                        }

                        double checkrsum = 0;
                        for (int i = 0; i < nbEle; ++i)
                                checkrsum+=data[i];
                        std::cout << "checkrsum: " << std::fixed  << checkrsum << std::endl;

                        double rsum = 0;
                        for (int i = 0; i < realNbEle; ++i)
                                rsum+=realData[i];
                        std::cout << "rsum: "  << std::fixed << rsum << std::endl;



                        std::cout << "writeDataSize: " << writeDataSize << std::endl;

                        std::cout << "procOffsets: ";
                        for(auto i = 0; i < procLength.size(); ++i)
                                std::cout << procOffsets[i] << " ";
                        std::cout << std::endl;


                        /*Print the first 10 data values to check the correctness.*/
                        // int i;
                        // printf("reconstructed data = ");
                        // for(i=0;i<20;i++)
                        //      printf("%f ", data[i]);
                        // printf("\n");
                        // free(data);
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
        // hid_t offDataset, offDatatype, offDataspace;
    // herr_t offStatus;
    // offDataset = H5Dopen2(file, "/level_0/data:offsets=0", H5P_DEFAULT);
    // offDatatype = H5Dget_type(offDataset);
        // offDataspace = H5Dget_space(offDataset);
        // int offLen = H5Sget_simple_extent_npoints(offDataspace);
        // // std::cout << "offLen: " << offLen << std::endl;
        // long long offset[offLen];
    // offStatus = H5Dread(offDataset, offDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, offset);
        // // std::cout << "offset: ";
        // // for (int i = 0; i < offLen; ++i)
        // //   std::cout <<  offset[i] << " ";
        // // std::cout << std::endl;
        // offStatus = H5Dclose(offDataset);*/

        status = H5Pclose(dcpl);
        status = H5Dclose(dset);
        status = H5Fclose(file);

        std::vector<int> spList;
        int sp;
        sprintf(meta_name, "../meta/sp_%d.txt", LEVEL);
        FILE* spfile = fopen(meta_name, "r+");
        if (spfile == NULL) {
                printf("Could not open spfile!\n");
                return 1;
        }
        while (fscanf(spfile, "%d", &sp) == 1)
                spList.push_back(sp);
        fclose(spfile);
        spList.push_back(-1);
        std::cout << "spList: ";
        for(auto it = spList.begin(); it != spList.end(); ++it) {
                std::cout << *it << " ";
        }
        std::cout << std::endl;

        double *testArr = new double[XSIZE*YSIZE*ZSIZE]{1.0};
        for (int i = 0; i < XSIZE*YSIZE*ZSIZE; ++i)
                testArr[i] = 1;

        for (int j = 0; j < boxLen_1; ++j) {
                for (int i = 0; i < 6; ++i){
                        if (i < 3)
                                box_1[j][i] = box_1[j][i]/2;
                        else
                                box_1[j][i] = (box_1[j][i]+1)/2 - 1;
                }
                // std::cout << " " << std::endl;
        }

        for (int i = 0; i < boxLen_1; ++i) {
                        for (size_t z=box_1[i][2]; z<=box_1[i][5]; ++z)
                                for (size_t y=box_1[i][1]; y<=box_1[i][4]; ++y)
                                        for (size_t x=box_1[i][0]; x<=box_1[i][3]; ++x) {
                                                testArr[x + y*XSIZE + z*XSIZE*YSIZE] = 0;
                                        }
        }

        long long tempCnt = 0;
        long long cnt = 0;
        long long offset;

        // /*1d decomp*/
        // for (int j = 0; j < boxLen_0; ++j) {
        //      for (size_t z=box_0[j][2]; z<=box_0[j][5]; ++z)
        //              for (size_t y=box_0[j][1]; y<=box_0[j][4]; ++y)
        //                      for (size_t x=box_0[j][0]; x<=box_0[j][3]; ++x) {
        //                              if (testArr[x + y*64 + z*64*64] != 0) {
        //                                      // std::cout <<x + y*64 + z*64*64 << std::endl;
        //                                      testArr[x + y*64 + z*64*64] = realData[cnt];
        //                                      cnt++;
        //                              }
        //                      }
        // }

        /*nast decomp*/
        for (int a = 0; a < boxLen_0; ++a) {
                        if (a==0 || spList[a]!= spList[a-1]) {
                                offset = procOffsets[spList[a]];
                                // std::cout << "tempCnt: " << tempCnt <<std::endl;
                                cnt+=tempCnt;
                                tempCnt = 0;
                        }
                for (int z = 0; z < (box_0[a][5]-box_0[a][2]+1)/BSIZE; ++z)
                        for (int y = 0; y < (box_0[a][4]-box_0[a][1]+1)/BSIZE; ++y)
                                for (int x = 0; x < (box_0[a][3]-box_0[a][0]+1)/BSIZE; ++x)
                                        for (int k = box_0[a][2]+z*BSIZE; k < box_0[a][2]+z*BSIZE+BSIZE; ++k)
                                                for (int j =box_0[a][1]+y*BSIZE; j <box_0[a][1]+y*BSIZE+BSIZE; ++j)
                                                        for (int i = box_0[a][0]+x*BSIZE; i <box_0[a][0]+x*BSIZE+BSIZE; ++i){
                                                                if(testArr[i + j*XSIZE + k*XSIZE*YSIZE] != 0) {
                                                                        testArr[i + j*XSIZE + k*XSIZE*YSIZE] = realData[offset+tempCnt];
                                                                        tempCnt++;
                                                                }
                                                        }
        }
        cnt+=tempCnt;
        double nastSum = 0;
        for (int i = 0; i < XSIZE * YSIZE*ZSIZE;++i)
                nastSum+=testArr[i];
        std::cout<< "nast sum: "<< nastSum << std::endl;
    // /*stack decomp*/


        // std::vector<int> addList;
        // int add;
        /*remove add*/
        // for (int i=0; i<nProcs; ++i) {
        //      sprintf(meta_name, "../meta/a_%d_%d.txt", LEVEL, i);
        //      FILE* afile = fopen(meta_name, "r+");
        //      if (afile == NULL) {
        //              printf("Could not open afile!\n");
        //              return 1;
        //      }
        //      fscanf(afile, "%d", &add);
        //      addList.push_back(add);
        //      fclose(afile);
        // }
        // std::cout << "addList: ";
        // for(auto it = addList.begin(); it != addList.end(); ++it) {
        //      std::cout << *it << " ";
        // }
        // std::cout << std::endl;

        // /*stack decomp*/
        // std::vector<int> bigList;
        // int big;
        // for (int i=0; i<nProcs; ++i) {
        //      sprintf(meta_name, "../meta/s_%d_%d.txt", LEVEL, i);
        //      FILE* bigfile = fopen(meta_name, "r+");
        //      if (bigfile == NULL) {
        //              printf("Could not open bigfile!\n");
        //              return 1;
        //      }
        //      fscanf(bigfile, "%d", &big);
        //      bigList.push_back(big);
        //      fclose(bigfile);
        // }
        // std::cout << "bigList: ";
        // for(auto it = bigList.begin(); it != bigList.end(); ++it) {
        //      std::cout << *it << " ";
        // }
        // std::cout << std::endl;

        // size_t unitBlkSize = BSIZE*BSIZE*BSIZE;
        // size_t bigX, big2X;
        // size_t b2Size = BSIZE*BSIZE;
        // size_t cc=0;
        // size_t xx, yy, zz, bb, ii, jj, kk;
        // long long sum = 0;
        // for (int a = 0; a < boxLen_0; ++a) {
        //      if (a==0 || spList[a]!= spList[a-1]) {
        //              offset = procOffsets[spList[a]];
        //              cnt = 0;
        //      }
        //      bigX = bigList[spList[a]];
        //      big2X = bigX*bigX;
        //      // std::cout << "a: " << a << std::endl;
        //      // std::cout << "bigX: " << bigX << std::endl;
        //      for (int z = 0; z < (box_0[a][5]-box_0[a][2]+1)/BSIZE; ++z){
        //              for (int y = 0; y < (box_0[a][4]-box_0[a][1]+1)/BSIZE; ++y){
        //                      for (int x = 0; x < (box_0[a][3]-box_0[a][0]+1)/BSIZE; ++x){
        //                              // todo bb = 0
        //                              for (int k = box_0[a][2]+z*BSIZE; k < box_0[a][2]+z*BSIZE+BSIZE; ++k){
        //                                      for (int j =box_0[a][1]+y*BSIZE; j <box_0[a][1]+y*BSIZE+BSIZE; ++j){
        //                                              for (int i = box_0[a][0]+x*BSIZE; i <box_0[a][0]+x*BSIZE+BSIZE; ++i){
        //                                                      if(testArr[i + j*XSIZE + k*XSIZE*YSIZE] != 0) {
        //                                                              cc = cnt/(unitBlkSize);
        //                                                              zz = cc/big2X;
        //                                                              yy = (cc - zz*big2X)/bigX;
        //                                                              xx = cc -zz*big2X - yy*bigX;
        //                                                              bb = cnt - unitBlkSize*cc;
        //                                                              kk = bb/b2Size;
        //                                                              jj = (bb - kk*b2Size)/BSIZE;
        //                                                              ii = bb - kk*b2Size - jj*BSIZE;
        //                                                              if (testArr[i + j*XSIZE + k*XSIZE*YSIZE] != 1)
        //                                                                  std::cout << "************" <<testArr[i + j*XSIZE + k*XSIZE*YSIZE] <<std::endl;
        //                                                              // std::cout <<xx << " "<< yy << " "<< zz << " "<< ii << " "<< jj << " "<< kk << " "<< bigX << " "<< big2X << " "<< BSIZE << " "<<std::endl;
        //                                                              // std::cout <<i << " "<< j << " "<< k <<std::endl;

        //                                                              // std::cout <<(xx*BSIZE+ii) + (yy*BSIZE+jj)*BSIZE*bigX + (zz*BSIZE+kk)*big2X*b2Size <<std::endl;
        //                                                              testArr[i + j*XSIZE + k*XSIZE*YSIZE] = realData[offset+(xx*BSIZE+ii) + (yy*BSIZE+jj)*BSIZE*bigX + (zz*BSIZE+kk)*big2X*b2Size];
        //                                                              // std::cout << testArr[i + j*XSIZE + k*XSIZE*YSIZE] <<std::endl;

        //                                                              sum+=testArr[i + j*XSIZE + k*XSIZE*YSIZE];
        //                                                              cnt++;
        //                                                              // std::cout <<cnt <<std::endl;
        //                                                      }
        //                                              }
        //                                      }
        //                              }
        //                      }
        //              }
        //      }
        //      // if (spList[a] != spList[a+1]){
        //      //      // std::cout <<"33333" <<std::endl;
        //      //      cnt += addList[spList[a]];
        //      // }
        // }
        // std::cout << "sum: " << sum << std::endl;

        /*bs decomp*/
        // for (int j = 0; j < boxLen_0; ++j)
        //      for (size_t z=box_0[j][2]; z<=box_0[j][5]; ++z)
        //              for (size_t y=box_0[j][1]; y<=box_0[j][4]; ++y)
        //                      for (size_t x=box_0[j][0]; x<=box_0[j][3]; ++x) {
        //                              testArr[x + y*XSIZE + z*XSIZE*YSIZE] = realData[cnt];
        //                              cnt++;
        //                      }


        std::cout << "cnt: " << cnt << std::endl;
        std::ofstream tout("test.raw", std::ios::binary);
        tout.write(reinterpret_cast<const char*>(testArr), XSIZE*YSIZE*ZSIZE * sizeof(double));
        tout.close();

        free(realData);
        return 0;
}

