#include <AMReX_VisMF.H>
#include <AMReX_AsyncOut.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FPC.H>
#include <AMReX_FabArrayUtility.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#include "hdf5.h"

#ifdef AMREX_USE_HDF5_ZFP
#include "H5Zzfp_lib.h"
#include "H5Zzfp_props.h"
#endif

#ifdef AMREX_USE_HDF5_SZ
#include "H5Z_SZ.h"
#endif

#ifdef AMREX_USE_HDF5_SZ3
#include "hdf5_sz3/include/H5Z_SZ3.hpp"
#endif

#include <fstream>
#include <iomanip>

namespace amrex {

#ifdef AMREX_USE_HDF5_ASYNC
hid_t es_id_g = 0;
#endif

static int CreateWriteHDF5AttrDouble(hid_t loc, const char *name, hsize_t n, const double *data)
{
    herr_t ret;
    hid_t attr, attr_space;
    hsize_t dims = n;

    attr_space = H5Screate_simple(1, &dims, NULL);

    attr = H5Acreate(loc, name, H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        printf("%s: Error with H5Acreate [%s]\n", __func__, name);
        return -1;
    }

    ret  = H5Awrite(attr, H5T_NATIVE_DOUBLE, (void*)data);
    if (ret < 0) {
        printf("%s: Error with H5Awrite [%s]\n", __func__, name);
        return -1;
    }
    H5Sclose(attr_space);
    H5Aclose(attr);
    return 1;
}

static int CreateWriteHDF5AttrInt(hid_t loc, const char *name, hsize_t n, const int *data)
{
    herr_t ret;
    hid_t attr, attr_space;
    hsize_t dims = n;

    attr_space = H5Screate_simple(1, &dims, NULL);

    attr = H5Acreate(loc, name, H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        printf("%s: Error with H5Acreate [%s]\n", __func__, name);
        return -1;
    }

    ret  = H5Awrite(attr, H5T_NATIVE_INT, (void*)data);
    if (ret < 0) {
        printf("%s: Error with H5Awrite [%s]\n", __func__, name);
        return -1;
    }
    H5Sclose(attr_space);
    H5Aclose(attr);
    return 1;
}

static int CreateWriteHDF5AttrString(hid_t loc, const char *name, const char* str)
{
    hid_t attr, atype, space;
    herr_t ret;

    BL_ASSERT(name);
    BL_ASSERT(str);

    space = H5Screate(H5S_SCALAR);
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(str)+1);
    H5Tset_strpad(atype,H5T_STR_NULLTERM);
    attr = H5Acreate(loc, name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        printf("%s: Error with H5Acreate [%s]\n", __func__, name);
        return -1;
    }

    ret = H5Awrite(attr, atype, str);
    if (ret < 0) {
        printf("%s: Error with H5Awrite[%s]\n", __func__, name);
        return -1;
    }

    H5Tclose(atype);
    H5Sclose(space);
    H5Aclose(attr);

    return 1;
}

#ifdef BL_USE_MPI
static void SetHDF5fapl(hid_t fapl, MPI_Comm comm)
#else
static void SetHDF5fapl(hid_t fapl)
#endif
{
#ifdef BL_USE_MPI
    H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL);

    // Alignment and metadata block size
    int alignment = 16 * 1024 * 1024;
    int blocksize =  4 * 1024 * 1024;
    H5Pset_alignment(fapl, alignment, alignment);
    H5Pset_meta_block_size(fapl, blocksize);

    // Collective metadata ops
    H5Pset_coll_metadata_write(fapl, true);
    H5Pset_all_coll_metadata_ops(fapl, true);

    // Defer cache flush
    H5AC_cache_config_t cache_config;
    cache_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    H5Pget_mdc_config(fapl, &cache_config);
    cache_config.set_initial_size = 1;
    cache_config.initial_size = 16 * 1024 * 1024;
    cache_config.evictions_enabled = 0;
    cache_config.incr_mode = H5C_incr__off;
    cache_config.flash_incr_mode = H5C_flash_incr__off;
    cache_config.decr_mode = H5C_decr__off;
    H5Pset_mdc_config (fapl, &cache_config);
#else
    H5Pset_fapl_sec2(fapl);
#endif

}

static void
WriteGenericPlotfileHeaderHDF5 (hid_t fid,
                               int nlevels,
                               const Vector<const MultiFab*>& mf,
                               const Vector<BoxArray> &bArray,
                               const Vector<std::string> &varnames,
                               const Vector<Geometry> &geom,
                               Real time,
                               const Vector<int> &level_steps,
                               const Vector<IntVect> &ref_ratio,
                               const std::string &versionName,
                               const std::string &levelPrefix,
                               const std::string &mfPrefix,
                               const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteGenericPlotfileHeaderHDF5()");

    BL_ASSERT(nlevels <= bArray.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());

    int finest_level(nlevels - 1);

    CreateWriteHDF5AttrString(fid, "version_name", versionName.c_str());
    CreateWriteHDF5AttrString(fid, "plotfile_type", "VanillaHDF5");

    int ncomp = varnames.size();
    CreateWriteHDF5AttrInt(fid, "num_components", 1, &ncomp);

    char comp_name[32];
    for (int ivar = 0; ivar < varnames.size(); ++ivar) {
        sprintf(comp_name, "component_%d", ivar);
        CreateWriteHDF5AttrString(fid, comp_name, varnames[ivar].c_str());
    }

    int ndim = AMREX_SPACEDIM;
    CreateWriteHDF5AttrInt(fid, "dim", 1, &ndim);
    double cur_time = (double)time;
    CreateWriteHDF5AttrDouble(fid, "time", 1, &cur_time);
    CreateWriteHDF5AttrInt(fid, "finest_level", 1, &finest_level);


    int coord = (int) geom[0].Coord();
    CreateWriteHDF5AttrInt(fid, "coordinate_system", 1, &coord);

    hid_t grp;
    char level_name[128];
    double lo[AMREX_SPACEDIM], hi[AMREX_SPACEDIM], cellsizes[AMREX_SPACEDIM];

    // For VisIt Chombo plot
    CreateWriteHDF5AttrInt(fid, "num_levels", 1, &nlevels);
    grp = H5Gcreate(fid, "Chombo_global", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    CreateWriteHDF5AttrInt(grp, "SpaceDim", 1, &ndim);
    H5Gclose(grp);

    hid_t comp_dtype;

    comp_dtype = H5Tcreate (H5T_COMPOUND, 2 * AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_j", 3 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_k", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 3 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_j", 4 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_k", 5 * sizeof(int), H5T_NATIVE_INT);
    }

    for (int level = 0; level <= finest_level; ++level) {
        sprintf(level_name, "level_%d", level);
        /* sprintf(level_name, "%s%d", levelPrefix.c_str(), level); */
        grp = H5Gcreate(fid, level_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (grp < 0) {
            std::cout << "H5Gcreate [" << level_name << "] failed!" << std::endl;
            continue;
        }

        int ratio = 1;
        if (ref_ratio.size() > 0)
            ratio = (level == finest_level)? 1: ref_ratio[level][0];

        CreateWriteHDF5AttrInt(grp, "ref_ratio", 1, &ratio);

        for (int k = 0; k < AMREX_SPACEDIM; ++k) {
            cellsizes[k] = (double)geom[level].CellSize()[k];
        }
        // Visit has issues with vec_dx, and is ok with a single "dx" value
        CreateWriteHDF5AttrDouble(grp, "Vec_dx", AMREX_SPACEDIM, cellsizes);
        // For VisIt Chombo plot
        CreateWriteHDF5AttrDouble(grp, "dx", 1, &cellsizes[0]);

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            lo[i] = (double)geom[level].ProbLo(i);
            hi[i] = (double)geom[level].ProbHi(i);
        }
        CreateWriteHDF5AttrDouble(grp, "prob_lo", AMREX_SPACEDIM, lo);
        CreateWriteHDF5AttrDouble(grp, "prob_hi", AMREX_SPACEDIM, hi);

        int domain[AMREX_SPACEDIM*2];
        Box tmp(geom[level].Domain());
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            domain[i] = tmp.smallEnd(i);
            domain[i+AMREX_SPACEDIM] = tmp.bigEnd(i);
        }

        hid_t aid = H5Screate(H5S_SCALAR);
        hid_t domain_attr = H5Acreate(grp, "prob_domain", comp_dtype, aid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(domain_attr, comp_dtype, domain);
        H5Aclose(domain_attr);
        H5Sclose(aid);

        int type[AMREX_SPACEDIM];
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            type[i] = (int)geom[level].Domain().ixType().test(i) ? 1 : 0;
        }
        CreateWriteHDF5AttrInt(grp, "domain_type", AMREX_SPACEDIM, type);

        CreateWriteHDF5AttrInt(grp, "steps", 1, &level_steps[level]);

        int ngrid = bArray[level].size();
        CreateWriteHDF5AttrInt(grp, "ngrid", 1, &ngrid);
        cur_time = (double)time;
        CreateWriteHDF5AttrDouble(grp, "time", 1, &cur_time);

        int ngrow = mf[level]->nGrow();
        CreateWriteHDF5AttrInt(grp, "ngrow", 1, &ngrow);

        /* hsize_t npts = ngrid*AMREX_SPACEDIM*2; */
        /* double *realboxes = new double [npts]; */
        /* for (int i = 0; i < bArray[level].size(); ++i) */
        /* { */
        /*     const Box &b(bArray[level][i]); */
        /*     RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo()); */
        /*     for (int n = 0; n < AMREX_SPACEDIM; ++n) { */
        /*         /1* HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n'; *1/ */
        /*         realboxes[i*AMREX_SPACEDIM*2 + n] = loc.lo(n); */
        /*         realboxes[i*AMREX_SPACEDIM*2 + AMREX_SPACEDIM + n] = loc.hi(n); */
        /*     } */
        /* } */
        /* CreateWriteDsetDouble(grp, "Boxes", npts, realboxes); */
        /* delete [] realboxes; */

        H5Gclose(grp);
    }

    H5Tclose(comp_dtype);
}

#ifdef AMREX_USE_HDF5_ASYNC
void async_vol_es_wait_close()
{
    size_t num_in_progress;
    hbool_t op_failed;
    if (es_id_g != 0) {
        H5ESwait(es_id_g, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
        if (num_in_progress != 0)
            std::cout << "After H5ESwait, still has async operations in progress!" << std::endl;
        H5ESclose(es_id_g);
        es_id_g = 0;
        /* std::cout << "es_id_g closed!" << std::endl; */
    }
    return;
}
static void async_vol_es_wait()
{
    size_t num_in_progress;
    hbool_t op_failed;
    if (es_id_g != 0) {
        H5ESwait(es_id_g, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
        if (num_in_progress != 0)
            std::cout << "After H5ESwait, still has async operations in progress!" << std::endl;
    }
    return;
}
#endif

void WriteMultiLevelPlotfileHDF5SingleDset (const std::string& plotfilename,
                                            int nlevels,
                                            const Vector<const MultiFab*>& mf,
                                            const Vector<std::string>& varnames,
                                            const Vector<Geometry>& geom,
                                            Real time,
                                            const Vector<int>& level_steps,
                                            const Vector<IntVect>& ref_ratio,
                                            const std::string &compression,
                                            const std::string &versionName,
                                            const std::string &levelPrefix,
                                            const std::string &mfPrefix,
                                            const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfileHDF5SingleDset");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int myProc(ParallelDescriptor::MyProc());
    int nProcs(ParallelDescriptor::NProcs());

#ifdef AMREX_USE_HDF5_ASYNC
    // For HDF5 async VOL, block and wait previous tasks have all completed
    if (es_id_g != 0) {
        async_vol_es_wait();
    }
    else {
        ExecOnFinalize(async_vol_es_wait_close);
        es_id_g = H5EScreate();
    }
#endif

    herr_t  ret;
    int finest_level = nlevels-1;
    int ncomp = mf[0]->nComp();
    /* double total_write_start_time(ParallelDescriptor::second()); */
    std::string filename(plotfilename + ".h5");

    // Write out root level metadata
    hid_t fapl, dxpl_col, dxpl_ind, dcpl_id, fid, grp;

    if(ParallelDescriptor::IOProcessor()) {
        BL_PROFILE_VAR("H5writeMetadata", h5dwm);
        // Create the HDF5 file
        fid = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (fid < 0)
            FileOpenFailed(filename.c_str());

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        WriteGenericPlotfileHeaderHDF5(fid, nlevels, mf, boxArrays, varnames, geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
        H5Fclose(fid);
        BL_PROFILE_VAR_STOP(h5dwm);
    }

    ParallelDescriptor::Barrier();

    hid_t babox_id;
    babox_id = H5Tcreate (H5T_COMPOUND, 2 * AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
        H5Tinsert (babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_i", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
        H5Tinsert (babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_i", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_j", 3 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
        H5Tinsert (babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "lo_k", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_i", 3 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_j", 4 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_k", 5 * sizeof(int), H5T_NATIVE_INT);
    }

    hid_t center_id = H5Tcreate (H5T_COMPOUND, AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
        H5Tinsert (center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
        H5Tinsert (center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (center_id, "j", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
        H5Tinsert (center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (center_id, "j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (center_id, "k", 2 * sizeof(int), H5T_NATIVE_INT);
    }

    fapl = H5Pcreate (H5P_FILE_ACCESS);
    dxpl_col = H5Pcreate(H5P_DATASET_XFER);
    dxpl_ind = H5Pcreate(H5P_DATASET_XFER);

#ifdef BL_USE_MPI
    SetHDF5fapl(fapl, ParallelDescriptor::Communicator());
    H5Pset_dxpl_mpio(dxpl_col, H5FD_MPIO_COLLECTIVE);
#else
    SetHDF5fapl(fapl);
#endif

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);

#if (defined AMREX_USE_HDF5_ZFP) || (defined AMREX_USE_HDF5_SZ) || (defined AMREX_USE_HDF5_SZ3)
    const char *chunk_env = NULL;
    std::string mode_env, value_env;
    double comp_value = -1.0;
    hsize_t chunk_dim[1] = {1024};

    chunk_env = getenv("HDF5_CHUNK_SIZE");
    if (chunk_env != NULL)
        chunk_dim[0] = atoi(chunk_env);

    H5Pset_chunk(dcpl_id, 1, chunk_dim);
    H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_INCR);

    std::string::size_type pos = compression.find('@');
    if (pos != std::string::npos) {
        mode_env = compression.substr(0, pos);
        value_env = compression.substr(pos+1);
        if (!value_env.empty()) {
            comp_value = atof(value_env.c_str());
        }
    }

#ifdef AMREX_USE_HDF5_ZFP
    pos = compression.find("ZFP");
    if (pos != std::string::npos) {
        ret = H5Z_zfp_initialize();
        if (ret < 0) amrex::Abort("ZFP initialize failed!");
    }
#endif

// #ifdef AMREX_USE_HDF5_SZ
//     pos = compression.find("SZ");
//     if (pos != std::string::npos) {
//         ret = H5Z_SZ_Init((char*)value_env.c_str());
//         if (ret < 0) {
//             std::cout << "SZ config file:" << value_env.c_str() << std::endl;
//             amrex::Abort("SZ initialize failed, check SZ config file!");
//         }
//     }
// #endif

    if (!mode_env.empty() && mode_env != "None") {
        if (mode_env == "ZLIB")
            H5Pset_deflate(dcpl_id, (int)comp_value);
#ifdef AMREX_USE_HDF5_ZFP
        else if (mode_env == "ZFP_RATE")
            H5Pset_zfp_rate(dcpl_id, comp_value);
        else if (mode_env == "ZFP_PRECISION")
            H5Pset_zfp_precision(dcpl_id, (unsigned int)comp_value);
        else if (mode_env == "ZFP_ACCURACY")
            H5Pset_zfp_accuracy(dcpl_id, comp_value);
        else if (mode_env == "ZFP_REVERSIBLE")
            H5Pset_zfp_reversible(dcpl_id);
        else if (mode_env == "ZLIB")
            H5Pset_deflate(dcpl_id, (int)comp_value);
#endif

        if (ParallelDescriptor::MyProc() == 0) {
            std::cout << "\nHDF5 plotfile using " << mode_env << std::endl;
            // << ", " <<
            //     value_env << ", " << chunk_dim[0] << std::endl;
        }
    }
#endif

    BL_PROFILE_VAR("H5writeAllLevel", h5dwd);

    // All process open the file
#ifdef AMREX_USE_HDF5_ASYNC
    // Only use async for writing actual data
    fid = H5Fopen_async(filename.c_str(), H5F_ACC_RDWR, fapl, es_id_g);
#else
    fid = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
#endif
    if (fid < 0)
        FileOpenFailed(filename.c_str());

    auto whichRD = FArrayBox::getDataDescriptor();
    bool doConvert(*whichRD != FPC::NativeRealDescriptor());
    int whichRDBytes(whichRD->numBytes());

    // Write data for each level
    char level_name[32];
    for (int level = 0; level <= finest_level; ++level) {
        hid_t  dcpl_id_lev;
        dcpl_id_lev = H5Pcopy(dcpl_id);
        sprintf(level_name, "level_%d", level);
#ifdef AMREX_USE_HDF5_ASYNC
        grp = H5Gopen_async(fid, level_name, H5P_DEFAULT, es_id_g);
#else
        grp = H5Gopen(fid, level_name, H5P_DEFAULT);
#endif
        if (grp < 0) { std::cout << "H5Gopen [" << level_name << "] failed!" << std::endl; break; }

        // Get the boxes assigned to all ranks and calculate their offsets and sizes
        Vector<int> procMap = mf[level]->DistributionMap().ProcessorMap();
        const BoxArray& grids = mf[level]->boxArray();
        hid_t boxdataset, boxdataspace;
        hid_t offsetdataset, offsetdataspace;
        hid_t centerdataset, centerdataspace;
        std::string bdsname("boxes");
        std::string odsname("data:offsets=0");
        std::string centername("boxcenter");
        std::string dataname("data:datatype=0");
        hsize_t  flatdims[1];
        flatdims[0] = grids.size();

        flatdims[0] = grids.size();
        boxdataspace = H5Screate_simple(1, flatdims, NULL);

#ifdef AMREX_USE_HDF5_ASYNC
        boxdataset = H5Dcreate_async(grp, bdsname.c_str(), babox_id, boxdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_id_g);
#else
        boxdataset = H5Dcreate(grp, bdsname.c_str(), babox_id, boxdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
        if (boxdataset < 0) { std::cout << "H5Dcreate [" << bdsname << "] failed!" << std::endl; break; }


        // Create a boxarray sorted by rank
        std::map<int, Vector<Box> > gridMap;
        for(int i(0); i < grids.size(); ++i) {
            int gridProc(procMap[i]);
            Vector<Box> &boxesAtProc = gridMap[gridProc];
            boxesAtProc.push_back(grids[i]);
        }
        BoxArray sortedGrids(grids.size());
        Vector<int> sortedProcs(grids.size());
        int bIndex(0);
        for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
            int proc = it->first;
            Vector<Box> &boxesAtProc = it->second;
            for(int ii(0); ii < boxesAtProc.size(); ++ii) {
                sortedGrids.set(bIndex, boxesAtProc[ii]);
                sortedProcs[bIndex] = proc;
                ++bIndex;
            }
        }

        // //dcdc-output procMap & sortedProcs
        // if(ParallelDescriptor::IOProcessor()) {
        //     std::cout << "procMap: ";
        //     for (auto i = procMap.begin(); i != procMap.end(); ++i)
        //         std::cout << *i << " ";
        //     std::cout << std::endl;

        //     std::cout << "sortedProcs: ";
        //     for (auto i = sortedProcs.begin(); i != sortedProcs.end(); ++i)
        //         std::cout << *i << " ";
        //     std::cout << std::endl;
        // }


        hsize_t  oflatdims[1];
        oflatdims[0] = sortedGrids.size() + 1;
        offsetdataspace = H5Screate_simple(1, oflatdims, NULL);
#ifdef AMREX_USE_HDF5_ASYNC
        offsetdataset = H5Dcreate_async(grp, odsname.c_str(), H5T_NATIVE_LLONG, offsetdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_id_g);
#else
        offsetdataset = H5Dcreate(grp, odsname.c_str(), H5T_NATIVE_LLONG, offsetdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
        if(offsetdataset < 0) { std::cout << "create offset dataset failed! ret = " << offsetdataset << std::endl; break;}

        hsize_t centerdims[1];
        centerdims[0]   = sortedGrids.size() ;
        centerdataspace = H5Screate_simple(1, centerdims, NULL);
#ifdef AMREX_USE_HDF5_ASYNC
        centerdataset = H5Dcreate_async(grp, centername.c_str(), center_id, centerdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_id_g);
#else
        centerdataset = H5Dcreate(grp, centername.c_str(), center_id, centerdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
        if(centerdataset < 0) { std::cout << "Create center dataset failed! ret = " << centerdataset << std::endl; break;}

        Vector<unsigned long long> offsets(sortedGrids.size() + 1, 0);
        unsigned long long currentOffset(0L);
        for(int b(0); b < sortedGrids.size(); ++b) {
            offsets[b] = currentOffset;
            currentOffset += sortedGrids[b].numPts() * ncomp;
        }
        offsets[sortedGrids.size()] = currentOffset;


        Vector<unsigned long long> procOffsets(nProcs, 0);
        Vector<unsigned long long> procBufferSize(nProcs, 0);
        Vector<unsigned long long> realProcBufferSize(nProcs, 0);
        unsigned long long totalOffset(0);
        for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
            int proc = it->first;
            Vector<Box> &boxesAtProc = it->second;
            procOffsets[proc] = totalOffset;
            procBufferSize[proc] = 0L;

            for(int b(0); b < boxesAtProc.size(); ++b) {
                procBufferSize[proc] += boxesAtProc[b].numPts() * ncomp;
            }


            totalOffset += procBufferSize[proc];
            /* if (level == 2) { */
            /*     fprintf(stderr, "Rank %d: level %d, proc %d, offset %ld, size %ld, all size %ld\n", */
            /*             myProc, level, proc, procOffsets[proc], procBufferSize[proc], totalOffset); */
            /* } */
        }

        // //dcdc-output procOffsets
        if(ParallelDescriptor::IOProcessor()) {
            std::cout << "procOffsets: ";
            for (auto i = procOffsets.begin(); i != procOffsets.end(); ++i)
                std::cout << *i << " ";
            std::cout << std::endl;

        }

        //dcdc output procBufferSize
        // if(ParallelDescriptor::IOProcessor()) {
        //     for (int i(0); i < nProcs; ++i) {
        //         std::cout << "buffer size on proc " << i << ": " << procBufferSize[i] << std::endl;
        //     }
        // }

        // //dcdc output box
        // std::cout << "box on proc " << myProc << std::endl;
        // for (int i(0); i < gridMap[myProc].size(); ++i) {
        //     for(int j(0); j < AMREX_SPACEDIM; ++j) {
        //         std::cout << gridMap[myProc][i].smallEnd(j) << " " << gridMap[myProc][i].bigEnd(j) << std::endl;
        //     }
        // }

        //dcdc find maxBuf
        // unsigned long long maxBuf = *max_element(procBufferSize.begin(), procBufferSize.end());
        // if(ParallelDescriptor::IOProcessor())
        //     std::cout << "maxBuf: " << maxBuf << std::endl;
        // H5Pset_chunk(dcpl_id, 1, &maxBuf);

        // //dcdc try fill
        // Vector<unsigned long long> fillProcOffsets(nProcs, 0);
        // Vector<unsigned long long> fillProcBufferSize(nProcs, 0);
        // unsigned long long totalOffset(0);
        // for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
        //     int proc = it->first;
        //     Vector<Box> &boxesAtProc = it->second;
        //     fillProcOffsets[proc] = totalOffset;
        //     fillProcBufferSize[proc] = 0L;

        //     //dcdc try filling
        //     fillProcBufferSize[proc] += maxBuf;

        //     totalOffset += fillProcBufferSize[proc];
        // }


        if(ParallelDescriptor::IOProcessor()) {
            int vbCount(0);
            Vector<int> vbox(sortedGrids.size() * 2 * AMREX_SPACEDIM);
            Vector<int> centering(sortedGrids.size() * AMREX_SPACEDIM);
            for(int b(0); b < sortedGrids.size(); ++b) {
                for(int i(0); i < AMREX_SPACEDIM; ++i) {
                    vbox[(vbCount * 2 * AMREX_SPACEDIM) + i] = sortedGrids[b].smallEnd(i);
                    vbox[(vbCount * 2 * AMREX_SPACEDIM) + i + AMREX_SPACEDIM] = sortedGrids[b].bigEnd(i);
                    centering[vbCount * AMREX_SPACEDIM + i] = sortedGrids[b].ixType().test(i) ? 1 : 0;
                }
                ++vbCount;
            }

            // Only proc zero needs to write out this information
#ifdef AMREX_USE_HDF5_ASYNC
            ret = H5Dwrite_async(offsetdataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, dxpl_ind, &(offsets[0]), es_id_g);
#else
            ret = H5Dwrite(offsetdataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, dxpl_ind, &(offsets[0]));
#endif
            if(ret < 0) { std::cout << "Write offset dataset failed! ret = " << ret << std::endl; }

#ifdef AMREX_USE_HDF5_ASYNC
            ret = H5Dwrite_async(centerdataset, center_id, H5S_ALL, H5S_ALL, dxpl_ind, &(centering[0]), es_id_g);
#else
            ret = H5Dwrite(centerdataset, center_id, H5S_ALL, H5S_ALL, dxpl_ind, &(centering[0]));
#endif
            if(ret < 0) { std::cout << "Write center dataset failed! ret = " << ret << std::endl; }

#ifdef AMREX_USE_HDF5_ASYNC
            ret = H5Dwrite_async(boxdataset, babox_id, H5S_ALL, H5S_ALL, dxpl_ind, &(vbox[0]), es_id_g);
#else
            ret = H5Dwrite(boxdataset, babox_id, H5S_ALL, H5S_ALL, dxpl_ind, &(vbox[0]));
#endif
            if(ret < 0) { std::cout << "Write box dataset failed! ret = " << ret << std::endl; }
        } // end IOProcessor

        hsize_t hs_procsize[1], hs_allprocsize[1], ch_offset[1];

        ch_offset[0]       = procOffsets[myProc];          // ---- offset on this proc
        hs_procsize[0]     = procBufferSize[myProc];       // ---- size of buffer on this proc
        std::cout << " " << std::endl;
        // std::cout << "size of buffer on proc " << myProc << ": " << realProcBufferSize[myProc] << std::endl;
        std::cout << "size of buffer on proc " << myProc << ": " << procBufferSize[myProc] << std::endl;
        hs_allprocsize[0]  = offsets[sortedGrids.size()];  // ---- size of buffer on all procs
        //dcdc change total buf size
        // hs_allprocsize[0]     = maxBuf * nProcs ;       // ---- size of buffer on all procs

        hid_t dataspace    = H5Screate_simple(1, hs_allprocsize, NULL);
        hid_t memdataspace = H5Screate_simple(1, hs_procsize, NULL);

        // hsize_t test_hs_procsize[3] = {32,32,32};
        // hsize_t test_hs_allprocsize[3] = {64,64,64};
        // hid_t dataspace    = H5Screate_simple(3, test_hs_procsize, NULL);
        // hid_t memdataspace = H5Screate_simple(3, test_hs_allprocsize, NULL);

        /* fprintf(stderr, "Rank %d: level %d, offset %ld, size %ld, all size %ld\n", myProc, level, ch_offset[0], hs_procsize[0], hs_allprocsize[0]); */

        if (hs_procsize[0] == 0)
            H5Sselect_none(dataspace);
        else
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, ch_offset, NULL, hs_procsize, NULL);

        auto preFileTime0 = amrex::second();
        Vector<Real> a_buffer(procBufferSize[myProc], -1.0);
        const MultiFab* data;
        std::unique_ptr<MultiFab> mf_tmp;
        if (mf[level]->nGrowVect() != 0) {
            mf_tmp = std::make_unique<MultiFab>(mf[level]->boxArray(),
                                                mf[level]->DistributionMap(),
                                                mf[level]->nComp(), 0, MFInfo(),
                                                mf[level]->Factory());
            MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
            data = mf_tmp.get();
        } else {
            data = mf[level];
        }

        Long writeDataItems(0), writeDataSize(0);
        for(MFIter mfi(*data); mfi.isValid(); ++mfi) {
            const FArrayBox &fab = (*data)[mfi];
            writeDataItems = fab.box().numPts() * (*data).nComp();
            // std::cout << "fab.box().numPts(): " << fab.box().numPts() << std::endl;
            // if(doConvert) {
            //     RealDescriptor::convertFromNativeFormat(static_cast<void *> (a_buffer.dataPtr()+writeDataSize),
            //                                             writeDataItems, fab.dataPtr(), *whichRD);
            // } else {    // ---- copy from the fab
                memcpy(static_cast<void *> (a_buffer.dataPtr()+writeDataSize),
                       fab.dataPtr(), writeDataItems * sizeof(double));
            // }
            writeDataSize += writeDataItems;
            // std::cout << "fab.dataPtr()[0]: " << fab.dataPtr()[0] << std::endl;
        }
        // std::cout << "a_buffer[0]: " << a_buffer[0] << std::endl;
        // std::cout << "writeDataSize on proc " << myProc << ": " << writeDataSize << std::endl;

        auto  preFileTime = amrex::second() - preFileTime0;
        const int IOProc2        = ParallelDescriptor::IOProcessorNumber();
        ParallelDescriptor::ReduceRealMax(preFileTime,IOProc2);
        amrex::Print() << "pre time = " << preFileTime << "  seconds" << "\n\n";

        /*out a_buffer to .bin*/
        // if (myProc == 0)
        // {
        //     std::ofstream fout("aa.bin", std::ios::binary);
        //     fout.write((char*)&a_buffer[0], procBufferSize[myProc] * sizeof(double));
        //     fout.close();
        // }

        BL_PROFILE_VAR("H5DwriteData", h5dwg);

        double eb = 0.001;
        std::ifstream inputFile("eb.txt");
        if (inputFile.is_open()) {
            inputFile >> eb;
            inputFile.close();
        } else {
            std::cerr << "Unable to open the eb file." << std::endl;
        }

#ifdef AMREX_USE_HDF5_SZ
        if (mode_env == "SZ") {
            size_t cd_nelmts;
            unsigned int* cd_values = NULL;
            unsigned filter_config;
            SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 1, 10000, eb, 0, 1024);
            H5Pset_filter(dcpl_id_lev, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        }
#endif

#ifdef AMREX_USE_HDF5_SZ3
        if (mode_env == "SZ") {
            size_t cd_nelmts;
            unsigned int* cd_values = NULL;
            unsigned filter_config;
            SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 1, 1000, eb, 0, 1024);
            H5Pset_filter(dcpl_id_lev, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        }
#endif

#ifdef AMREX_USE_HDF5_ASYNC
        hid_t dataset = H5Dcreate_async(grp, dataname.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, dcpl_id_lev, H5P_DEFAULT, es_id_g);
#else
        hid_t dataset = H5Dcreate(grp, dataname.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, dcpl_id_lev, H5P_DEFAULT);
#endif
        if(dataset < 0)
            std::cout << ParallelDescriptor::MyProc() << "create data failed!  ret = " << dataset << std::endl;

        auto dPlotFileTime0 = amrex::second();

#ifdef AMREX_USE_HDF5_ASYNC
        ret = H5Dwrite_async(dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace, dxpl_col, a_buffer.dataPtr(), es_id_g);
#else
        ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace, dxpl_col, a_buffer.dataPtr());
#endif
        if(ret < 0) { std::cout << ParallelDescriptor::MyProc() << "Write data failed!  ret = " << ret << std::endl; break; }

        auto  dPlotFileTime = amrex::second() - dPlotFileTime0;
        const int IOProc        = ParallelDescriptor::IOProcessorNumber();
        ParallelDescriptor::ReduceRealMax(dPlotFileTime,IOProc);
        amrex::Print() << "real write time = " << dPlotFileTime << "  seconds" << "\n\n";

        BL_PROFILE_VAR_STOP(h5dwg);
        H5Pclose(dcpl_id_lev);
        H5Sclose(memdataspace);
        H5Sclose(dataspace);
        H5Sclose(offsetdataspace);
        H5Sclose(centerdataspace);
        H5Sclose(boxdataspace);

#ifdef AMREX_USE_HDF5_ASYNC
        H5Dclose_async(dataset, es_id_g);
        H5Dclose_async(offsetdataset, es_id_g);
        H5Dclose_async(centerdataset, es_id_g);
        H5Dclose_async(boxdataset, es_id_g);
        H5Gclose_async(grp, es_id_g);
#else
        H5Dclose(dataset);
        H5Dclose(offsetdataset);
        H5Dclose(centerdataset);
        H5Dclose(boxdataset);
        H5Gclose(grp);
#endif
    } // For group

    BL_PROFILE_VAR_STOP(h5dwd);

    H5Tclose(center_id);
    H5Tclose(babox_id);
    H5Pclose(fapl);
    H5Pclose(dxpl_col);
    H5Pclose(dxpl_ind);
    H5Pclose(dcpl_id);

#ifdef AMREX_USE_HDF5_ASYNC
    H5Fclose_async(fid, es_id_g);
#else
    H5Fclose(fid);
#endif

// //dcdc-start
//     for (int lev=finest_level; lev>=0; --lev)
//     {
//         MultiFab tempmf(mf[lev]->boxArray(),mf[lev]->DistributionMap(),mf[lev]->nComp(),0);
//         Real bogus_flag=-1e200;
//         if (lev < finest_level)
//         {
//             const BoxArray baf = BoxArray(tempmf.boxArray()).coarsen(ref_ratio[lev]);
//             for (MFIter mfi(tempmf); mfi.isValid(); ++mfi)
//             {
//                 FArrayBox& myFab = tempmf[mfi];
//                     std::vector< std::pair<int,Box> > isects = baf.intersections(tempmf.boxArray()[mfi.index()]);

//                 for (int ii = 0; ii < isects.size(); ii++)
//                     myFab.setVal(bogus_flag,isects[ii].second,0,ncomp);
//             }
//         }

//         long int idx=0;
//         // std::cout<<"Probhi: "<<geom[lev].ProbHi(0) <<  " " << geom[lev].ProbHi(1) << " " << geom[lev].ProbHi(2) <<std::endl;
//         double max[3] = { -FLT_MAX,-FLT_MAX,-FLT_MAX};
//         double min[3] = { FLT_MAX,FLT_MAX,FLT_MAX};
//         for (MFIter mfi(tempmf); mfi.isValid(); ++mfi)
//         {
//             const FArrayBox& fab = tempmf[mfi];
//              Array4<Real> const& fab_array = tempmf.array(mfi);
//              int ncomp = tempmf.nComp();
//             const Box& box = mfi.validbox();
//              FArrayBox fab2;
//              fab2.resize(box,ncomp);

//             const Dim3 lo = amrex::lbound(box);
//             const Dim3 hi = amrex::ubound(box);

//             const auto dx = geom[lev].CellSize();
//             RealVect boxSizeH(geom[lev].ProbHi());
//             // if(verbose)
//             // {
//                 // std::cout<<"level"<<lev<<"dx: "<<dx[0]<<"dy: "<<dx[1]<<dx[2]<<std::endl;
//                 // std::cout<<"Probhi"<<boxSizeH<<std::endl;
//             // }

//          const auto len = amrex::length(box);  // length of box

//          for (int z = lo.z; z <= hi.z; ++z) {
//              for (int y = lo.y; y <= hi.y; ++y) {
//                      for (int x = lo.x; x <= hi.x; ++x) {

//                 //idx assumes ncomp=1, is local to this fab
//                 //                         int idx = x+y*len.x+z*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
//                     if(fab_array(x,y,z,0) != bogus_flag) {

//                         Real tmpx=x*dx[0];
//                         Real tmpy=y*dx[1];
//                         Real tmpz=z*dx[2];
//                         Real tmpvol=dx[2]*dx[0]*dx[1];

//                         if (x > max[0])
//                                      max[0] = x;
//                         if (y > max[1])
//                                      max[1] = y;
//                         if (z > max[2])
//                                      max[2] = z;
//                         if (x < min[0])
//                                      min[0] = x;
//                         if (y < min[1])
//                                      min[1] = y;
//                         if (z < min[2])
//                                      min[2] = z;


//                         // for (int n = 0; n < comps.size(); ++n) {
//                         //     ofs_mpi.write((char *) &(fab_array(x,y,z,comps[n])), sizeof(amrex::Real));
//                         // }

//                     }
//                      }
//              }
//          }
//         std::cout<< "box on proc " << myProc << std::endl;
//         std::cout<< "low: " <<  lo.x<< " " << lo.y << " " << lo.z<< std::endl;
//         std::cout<< "high: " <<  hi.x<< " " << hi.y << " " << hi.z<< std::endl;
//         // std::cout<< "range: " <<  len.x<< " " << len.y << " " << len.z<< std::endl;
//         }
//         // std::cout<< "pro low: " <<  min[0]<< " " << min[1] << " " << min[2]<< std::endl;
//         // std::cout<< "pro high: " <<  max[0]<< " " << max[1] << " " << max[2]<< std::endl;
//         // std::cout<< "pro range: " <<  max[0]-min[0] << " " << max[1]-min[1] << " " << max[2]-min[2] << std::endl;
//         // std::cout<< " " << std::endl;

//     }

} // WriteMultiLevelPlotfileHDF5SingleDset

void WriteMultiLevelPlotfileHDF5MultiDset (const std::string& plotfilename,
                                           int nlevels,
                                           const Vector<const MultiFab*>& mf,
                                           const Vector<std::string>& varnames,
                                           const Vector<Geometry>& geom,
                                           Real time,
                                           const Vector<int>& level_steps,
                                           const Vector<IntVect>& ref_ratio,
                                           const std::string &compression,
                                           const std::string &versionName,
                                           const std::string &levelPrefix,
                                           const std::string &mfPrefix,
                                           const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfileHDF5MultiDset");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int myProc(ParallelDescriptor::MyProc());
    int nProcs(ParallelDescriptor::NProcs());

#ifdef AMREX_USE_HDF5_ASYNC
    // For HDF5 async VOL, block and wait previous tasks have all completed
    if (es_id_g != 0) {
        async_vol_es_wait();
    }
    else {
        ExecOnFinalize(async_vol_es_wait_close);
        es_id_g = H5EScreate();
    }
#endif

    herr_t  ret;
    int finest_level = nlevels-1;
    int ncomp = mf[0]->nComp();
    /* double total_write_start_time(ParallelDescriptor::second()); */
    std::string filename(plotfilename + ".h5");

    // Write out root level metadata
    hid_t fapl, dxpl_col, dxpl_ind, fid, grp, dcpl_id, dcpl_id_lev;

    if(ParallelDescriptor::IOProcessor()) {
        BL_PROFILE_VAR("H5writeMetadata", h5dwm);
        // Create the HDF5 file
        fid = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (fid < 0)
            FileOpenFailed(filename.c_str());

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        WriteGenericPlotfileHeaderHDF5(fid, nlevels, mf, boxArrays, varnames, geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
        H5Fclose(fid);
        BL_PROFILE_VAR_STOP(h5dwm);
    }

    ParallelDescriptor::Barrier();

    hid_t babox_id;
    babox_id = H5Tcreate (H5T_COMPOUND, 2 * AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
        H5Tinsert (babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_i", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
        H5Tinsert (babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_i", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_j", 3 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
        H5Tinsert (babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "lo_k", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_i", 3 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_j", 4 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (babox_id, "hi_k", 5 * sizeof(int), H5T_NATIVE_INT);
    }

    hid_t center_id = H5Tcreate (H5T_COMPOUND, AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
        H5Tinsert (center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
        H5Tinsert (center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (center_id, "j", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
        H5Tinsert (center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (center_id, "j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (center_id, "k", 2 * sizeof(int), H5T_NATIVE_INT);
    }

    fapl = H5Pcreate (H5P_FILE_ACCESS);
    dxpl_col = H5Pcreate(H5P_DATASET_XFER);
    dxpl_ind = H5Pcreate(H5P_DATASET_XFER);

#ifdef BL_USE_MPI
    SetHDF5fapl(fapl, ParallelDescriptor::Communicator());
    H5Pset_dxpl_mpio(dxpl_col, H5FD_MPIO_COLLECTIVE);
#else
    SetHDF5fapl(fapl);
#endif

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);

#if (defined AMREX_USE_HDF5_ZFP) || (defined AMREX_USE_HDF5_SZ) || (defined AMREX_USE_HDF5_SZ3)
    const char *chunk_env = NULL;
    std::string mode_env, value_env;
    double comp_value = -1.0;
    hsize_t chunk_dim = 32768;

    chunk_env = getenv("HDF5_CHUNK_SIZE");
    if (chunk_env != NULL)
        chunk_dim = atoi(chunk_env);

    H5Pset_chunk(dcpl_id, 1, &chunk_dim);
    H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_INCR);

    std::string::size_type pos = compression.find('@');
    if (pos != std::string::npos) {
        mode_env = compression.substr(0, pos);
        value_env = compression.substr(pos+1);
        if (!value_env.empty()) {
            comp_value = atof(value_env.c_str());
        }
    }

#ifdef AMREX_USE_HDF5_ZFP
    pos = compression.find("ZFP");
    if (pos != std::string::npos) {
        ret = H5Z_zfp_initialize();
        if (ret < 0) amrex::Abort("ZFP initialize failed!");
    }
#endif

// #ifdef AMREX_USE_HDF5_SZ
//     pos = compression.find("SZ");
//     if (pos != std::string::npos) {
//         ret = H5Z_SZ_Init((char*)value_env.c_str());
//         if (ret < 0) amrex::Abort("ZFP initialize failed, check SZ config file!");
//     }
// #endif

    if (!mode_env.empty() && mode_env != "None") {
        if (mode_env == "ZLIB")
            H5Pset_deflate(dcpl_id, (int)comp_value);
#ifdef AMREX_USE_HDF5_ZFP
        else if (mode_env == "ZFP_RATE")
            H5Pset_zfp_rate(dcpl_id, comp_value);
        else if (mode_env == "ZFP_PRECISION")
            H5Pset_zfp_precision(dcpl_id, (unsigned int)comp_value);
        else if (mode_env == "ZFP_ACCURACY")
            H5Pset_zfp_accuracy(dcpl_id, comp_value);
        else if (mode_env == "ZFP_REVERSIBLE")
            H5Pset_zfp_reversible(dcpl_id);
#endif

        if (ParallelDescriptor::MyProc() == 0) {
            std::cout << "\nHDF5 checkpoint using " << mode_env << ", " <<
                value_env << ", " << chunk_dim << std::endl;
        }
    }
#endif

    BL_PROFILE_VAR("H5writeAllLevel", h5dwd);

    // All process open the file
#ifdef AMREX_USE_HDF5_ASYNC
    // Only use async for writing actual data
    fid = H5Fopen_async(filename.c_str(), H5F_ACC_RDWR, fapl, es_id_g);
#else
    fid = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
#endif
    if (fid < 0)
        FileOpenFailed(filename.c_str());

    auto whichRD = FArrayBox::getDataDescriptor();
    bool doConvert(*whichRD != FPC::NativeRealDescriptor());
    int whichRDBytes(whichRD->numBytes());

    // Write data for each level
    char level_name[32];

    for (int level = 0; level <= finest_level; ++level) {
        sprintf(level_name, "level_%d", level);
#ifdef AMREX_USE_HDF5_ASYNC
        grp = H5Gopen_async(fid, level_name, H5P_DEFAULT, es_id_g);
#else
        grp = H5Gopen(fid, level_name, H5P_DEFAULT);
#endif
        if (grp < 0) { std::cout << "H5Gopen [" << level_name << "] failed!" << std::endl; break; }

        // Get the boxes assigned to all ranks and calculate their offsets and sizes
        Vector<int> procMap = mf[level]->DistributionMap().ProcessorMap();
        const BoxArray& grids = mf[level]->boxArray();
        hid_t boxdataset, boxdataspace;
        hid_t offsetdataset, offsetdataspace;
        hid_t centerdataset, centerdataspace;
        std::string bdsname("boxes");
        std::string odsname("data:offsets=0");
        std::string centername("boxcenter");
        hsize_t  flatdims[1];
        flatdims[0] = grids.size();

        flatdims[0] = grids.size();
        boxdataspace = H5Screate_simple(1, flatdims, NULL);

#ifdef AMREX_USE_HDF5_ASYNC
        boxdataset = H5Dcreate_async(grp, bdsname.c_str(), babox_id, boxdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_id_g);
#else
        boxdataset = H5Dcreate(grp, bdsname.c_str(), babox_id, boxdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
        if (boxdataset < 0) { std::cout << "H5Dcreate [" << bdsname << "] failed!" << std::endl; break; }

        // Create a boxarray sorted by rank
        std::map<int, Vector<Box> > gridMap;
        for(int i(0); i < grids.size(); ++i) {
            int gridProc(procMap[i]);
            Vector<Box> &boxesAtProc = gridMap[gridProc];
            boxesAtProc.push_back(grids[i]);
        }
        BoxArray sortedGrids(grids.size());
        Vector<int> sortedProcs(grids.size());
        int bIndex(0);
        for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
            int proc = it->first;
            Vector<Box> &boxesAtProc = it->second;
            for(int ii(0); ii < boxesAtProc.size(); ++ii) {
                sortedGrids.set(bIndex, boxesAtProc[ii]);
                sortedProcs[bIndex] = proc;
                ++bIndex;
            }
        }

        hsize_t  oflatdims[1];
        oflatdims[0] = sortedGrids.size() + 1;
        offsetdataspace = H5Screate_simple(1, oflatdims, NULL);
#ifdef AMREX_USE_HDF5_ASYNC
        offsetdataset = H5Dcreate_async(grp, odsname.c_str(), H5T_NATIVE_LLONG, offsetdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_id_g);
#else
        offsetdataset = H5Dcreate(grp, odsname.c_str(), H5T_NATIVE_LLONG, offsetdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
        if(offsetdataset < 0) { std::cout << "create offset dataset failed! ret = " << offsetdataset << std::endl; break;}

        hsize_t centerdims[1];
        centerdims[0]   = sortedGrids.size() ;
        centerdataspace = H5Screate_simple(1, centerdims, NULL);
#ifdef AMREX_USE_HDF5_ASYNC
        centerdataset = H5Dcreate_async(grp, centername.c_str(), center_id, centerdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_id_g);
#else
        centerdataset = H5Dcreate(grp, centername.c_str(), center_id, centerdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
        if(centerdataset < 0) { std::cout << "Create center dataset failed! ret = " << centerdataset << std::endl; break;}

        Vector<unsigned long long> offsets(sortedGrids.size() + 1);
        unsigned long long currentOffset(0L);
        for(int b(0); b < sortedGrids.size(); ++b) {
            offsets[b] = currentOffset;
            /* currentOffset += sortedGrids[b].numPts() * ncomp; */
            currentOffset += sortedGrids[b].numPts();
        }
        offsets[sortedGrids.size()] = currentOffset;

        Vector<unsigned long long> procOffsets(nProcs);
        Vector<unsigned long long> procBufferSize(nProcs);
        unsigned long long totalOffset(0);
        for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
            int proc = it->first;
            Vector<Box> &boxesAtProc = it->second;
            procOffsets[proc] = totalOffset;
            procBufferSize[proc] = 0L;
            for(int b(0); b < boxesAtProc.size(); ++b) {
                /* procBufferSize[proc] += boxesAtProc[b].numPts() * ncomp; */
                procBufferSize[proc] += boxesAtProc[b].numPts();
            }
            totalOffset += procBufferSize[proc];
        }

        unsigned long long maxBuf = *max_element(procBufferSize.begin(), procBufferSize.end());
        if(ParallelDescriptor::IOProcessor())
            std::cout << "maxBuf: " << maxBuf << std::endl;
        maxBuf = 32768;
        H5Pset_chunk(dcpl_id, 1, &maxBuf);

        if(ParallelDescriptor::IOProcessor()) {
            int vbCount(0);
            Vector<int> vbox(sortedGrids.size() * 2 * AMREX_SPACEDIM);
            Vector<int> centering(sortedGrids.size() * AMREX_SPACEDIM);
            for(int b(0); b < sortedGrids.size(); ++b) {
                for(int i(0); i < AMREX_SPACEDIM; ++i) {
                    vbox[(vbCount * 2 * AMREX_SPACEDIM) + i] = sortedGrids[b].smallEnd(i);
                    vbox[(vbCount * 2 * AMREX_SPACEDIM) + i + AMREX_SPACEDIM] = sortedGrids[b].bigEnd(i);
                    centering[vbCount * AMREX_SPACEDIM + i] = sortedGrids[b].ixType().test(i) ? 1 : 0;
                }
                ++vbCount;
            }

            // Only proc zero needs to write out this information
#ifdef AMREX_USE_HDF5_ASYNC
            ret = H5Dwrite_async(offsetdataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, dxpl_ind, &(offsets[0]), es_id_g);
#else
            ret = H5Dwrite(offsetdataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, dxpl_ind, &(offsets[0]));
#endif
            if(ret < 0) { std::cout << "Write offset dataset failed! ret = " << ret << std::endl; }

#ifdef AMREX_USE_HDF5_ASYNC
            ret = H5Dwrite_async(centerdataset, center_id, H5S_ALL, H5S_ALL, dxpl_ind, &(centering[0]), es_id_g);
#else
            ret = H5Dwrite(centerdataset, center_id, H5S_ALL, H5S_ALL, dxpl_ind, &(centering[0]));
#endif
            if(ret < 0) { std::cout << "Write center dataset failed! ret = " << ret << std::endl; }

#ifdef AMREX_USE_HDF5_ASYNC
            ret = H5Dwrite_async(boxdataset, babox_id, H5S_ALL, H5S_ALL, dxpl_ind, &(vbox[0]), es_id_g);
#else
            ret = H5Dwrite(boxdataset, babox_id, H5S_ALL, H5S_ALL, dxpl_ind, &(vbox[0]));
#endif
            if(ret < 0) { std::cout << "Write box dataset failed! ret = " << ret << std::endl; }
        } // end IOProcessor

        hsize_t hs_procsize[1], hs_allprocsize[1], ch_offset[1];

        ch_offset[0]       = procOffsets[myProc];          // ---- offset on this proc
        hs_procsize[0]     = procBufferSize[myProc];       // ---- size of buffer on this proc
        std::cout << "size of buffer on proc " << myProc << ": " << hs_procsize[0] << std::endl;
        hs_allprocsize[0]  = offsets[sortedGrids.size()];  // ---- size of buffer on all procs

        hid_t dataspace    = H5Screate_simple(1, hs_allprocsize, NULL);
        hid_t memdataspace = H5Screate_simple(1, hs_procsize, NULL);

        if (hs_procsize[0] == 0)
            H5Sselect_none(dataspace);
        else
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, ch_offset, NULL, hs_procsize, NULL);

        Vector<Real> a_buffer(procBufferSize[myProc]*ncomp, -1.0);
        Vector<Real> a_buffer_ind(procBufferSize[myProc], -1.0);
        const MultiFab* data;
        std::unique_ptr<MultiFab> mf_tmp;
        if (mf[level]->nGrowVect() != 0) {
            mf_tmp = std::make_unique<MultiFab>(mf[level]->boxArray(),
                                                mf[level]->DistributionMap(),
                                                mf[level]->nComp(), 0, MFInfo(),
                                                mf[level]->Factory());
            MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
            data = mf_tmp.get();
        } else {
            data = mf[level];
        }

        hid_t dataset;
        char dataname[64];

        dcpl_id_lev = H5Pcopy(dcpl_id);
#ifdef AMREX_USE_HDF5_SZ
        if (mode_env == "SZ") {
            size_t cd_nelmts;
            unsigned int* cd_values = NULL;
            unsigned filter_config;
            // SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 0, 10, 0, 0, 0);
            SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 1, 0, 0.001, 0, hs_procsize[0]);
            H5Pset_filter(dcpl_id_lev, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        }
#endif

#ifdef AMREX_USE_HDF5_SZ3
        if (mode_env == "SZ") {
            size_t cd_nelmts;
            unsigned int* cd_values = NULL;
            unsigned filter_config;
            SZ_errConfigToCdArray(&cd_nelmts, &cd_values, 1, 0, 0.001, 0, hs_procsize[0]);
            H5Pset_filter(dcpl_id_lev, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        }
#endif

        BL_PROFILE_VAR("H5DwriteData", h5dwg);

        for (int jj = 0; jj < ncomp; jj++) {

            Long writeDataItems(0), writeDataSize(0);
            for(MFIter mfi(*data); mfi.isValid(); ++mfi) {
                const FArrayBox &fab = (*data)[mfi];
                writeDataItems = fab.box().numPts();
                if(doConvert) {
                    RealDescriptor::convertFromNativeFormat(static_cast<void *> (a_buffer.dataPtr()),
                                                            writeDataItems * ncomp, fab.dataPtr(), *whichRD);

                } else {    // ---- copy from the fab
                    memcpy(static_cast<void *> (a_buffer.dataPtr()),
                           fab.dataPtr(), writeDataItems * ncomp * whichRDBytes);
                }

                // Extract individual variable data
                memcpy(static_cast<void *> (a_buffer_ind.dataPtr() + writeDataSize),
                       static_cast<void *> (a_buffer.dataPtr() + jj*writeDataItems),
                       writeDataItems * whichRDBytes);

                writeDataSize += writeDataItems;
            }

            sprintf(dataname, "data:datatype=%d", jj);
#ifdef AMREX_USE_HDF5_ASYNC
            dataset = H5Dcreate_async(grp, dataname, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, dcpl_id_lev, H5P_DEFAULT, es_id_g);
            if(dataset < 0) std::cout << ParallelDescriptor::MyProc() << "create data failed!  ret = " << dataset << std::endl;
            ret = H5Dwrite_async(dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace, dxpl_col, a_buffer_ind.dataPtr(), es_id_g);
            if(ret < 0) { std::cout << ParallelDescriptor::MyProc() << "Write data failed!  ret = " << ret << std::endl; break; }
            H5Dclose_async(dataset, es_id_g);
#else
            dataset = H5Dcreate(grp, dataname, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, dcpl_id_lev, H5P_DEFAULT);
            if(dataset < 0) std::cout << ParallelDescriptor::MyProc() << "create data failed!  ret = " << dataset << std::endl;
            ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace, dxpl_col, a_buffer_ind.dataPtr());
            if(ret < 0) { std::cout << ParallelDescriptor::MyProc() << "Write data failed!  ret = " << ret << std::endl; break; }
            H5Dclose(dataset);
#endif
        }

        BL_PROFILE_VAR_STOP(h5dwg);
        H5Pclose(dcpl_id_lev);
        H5Sclose(memdataspace);
        H5Sclose(dataspace);
        H5Sclose(offsetdataspace);
        H5Sclose(centerdataspace);
        H5Sclose(boxdataspace);

#ifdef AMREX_USE_HDF5_ASYNC
        H5Dclose_async(offsetdataset, es_id_g);
        H5Dclose_async(centerdataset, es_id_g);
        H5Dclose_async(boxdataset, es_id_g);
        H5Gclose_async(grp, es_id_g);
#else
        H5Dclose(offsetdataset);
        H5Dclose(centerdataset);
        H5Dclose(boxdataset);
        H5Gclose(grp);
#endif
    } // For group

    BL_PROFILE_VAR_STOP(h5dwd);

    H5Tclose(center_id);
    H5Tclose(babox_id);
    H5Pclose(fapl);
    H5Pclose(dcpl_id);
    H5Pclose(dxpl_col);
    H5Pclose(dxpl_ind);

#ifdef AMREX_USE_HDF5_ASYNC
    H5Fclose_async(fid, es_id_g);
#else
    H5Fclose(fid);
#endif


} // WriteMultiLevelPlotfileHDF5MultiDset

void
WriteSingleLevelPlotfileHDF5 (const std::string& plotfilename,
                              const MultiFab& mf, const Vector<std::string>& varnames,
                              const Geometry& geom, Real time, int level_step,
                              const std::string &compression,
                              const std::string &versionName,
                              const std::string &levelPrefix,
                              const std::string &mfPrefix,
                              const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    WriteMultiLevelPlotfileHDF5(plotfilename, 1, mfarr, varnames, geomarr, time, level_steps, ref_ratio,
                                compression, versionName, levelPrefix, mfPrefix, extra_dirs);
}

void
WriteSingleLevelPlotfileHDF5SingleDset (const std::string& plotfilename,
                                        const MultiFab& mf, const Vector<std::string>& varnames,
                                        const Geometry& geom, Real time, int level_step,
                                        const std::string &compression,
                                        const std::string &versionName,
                                        const std::string &levelPrefix,
                                        const std::string &mfPrefix,
                                        const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    WriteMultiLevelPlotfileHDF5SingleDset(plotfilename, 1, mfarr, varnames, geomarr, time, level_steps, ref_ratio,
                                          compression, versionName, levelPrefix, mfPrefix, extra_dirs);
}

void
WriteSingleLevelPlotfileHDF5MultiDset (const std::string& plotfilename,
                                       const MultiFab& mf, const Vector<std::string>& varnames,
                                       const Geometry& geom, Real time, int level_step,
                                       const std::string &compression,
                                       const std::string &versionName,
                                       const std::string &levelPrefix,
                                       const std::string &mfPrefix,
                                       const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    WriteMultiLevelPlotfileHDF5MultiDset(plotfilename, 1, mfarr, varnames, geomarr, time, level_steps, ref_ratio,
                                         compression, versionName, levelPrefix, mfPrefix, extra_dirs);
}

void
WriteMultiLevelPlotfileHDF5 (const std::string &plotfilename,
                             int nlevels,
                             const Vector<const MultiFab*> &mf,
                             const Vector<std::string> &varnames,
                             const Vector<Geometry> &geom,
                             Real time,
                             const Vector<int> &level_steps,
                             const Vector<IntVect> &ref_ratio,
                             const std::string &compression,
                             const std::string &versionName,
                             const std::string &levelPrefix,
                             const std::string &mfPrefix,
                             const Vector<std::string>& extra_dirs)
{

    WriteMultiLevelPlotfileHDF5SingleDset(plotfilename, nlevels, mf, varnames, geom, time, level_steps, ref_ratio,
                                          compression, versionName, levelPrefix, mfPrefix, extra_dirs);
}

} // namespace amrex

