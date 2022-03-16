*     Developed by Maxwell Xu CAI (NAOC/KIAA)
*     For bug report/feedback, please email maxwellemail@gmail.com
*     The binary format is developed by Long Wang (May 11, 2015)
*
*     Version 2.0 -- Bug fixes & finetuning, June 12, 2014
*     Version 1.4 -- June 11, 2014
*     Version 1.0 -- Initial release, Nov 2013
*

      SUBROUTINE HDF5_INIT(h5_fname)
*
*
*       Initialize a new HDF5 file given by its file name.
*       --------------------------
*
*
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
*       SAVE Block
      INTEGER(HID_T) :: h5_file_id, h5_step, h5_group_id
      INTEGER(HID_T), DIMENSION(1) :: h5_dset_ids(256)
      INTEGER(HID_T) :: h5_vec_len
      REAL*8 :: h5_current_time
      CHARACTER(LEN=64) :: h5_file_name
      CHARACTER(LEN=64) :: h5_fname
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER*4 :: ERROR
      
*       Close any previously opened file (if any)
      IF (h5_file_id .GT. 0) CALL HDF5_CLOSE

*     Initialize FORTRAN interface.
      CALL h5open_f(ERROR) 

*     Create a new file using default properties.
      h5_file_name = h5_fname
      CALL h5fcreate_f(TRIM(h5_file_name), H5F_ACC_TRUNC_F, h5_file_id, 
     &      ERROR)

      h5_step = -1
      print*, 'HDF5 Init ID:',h5_file_id,' Snapshot file: ',h5_file_name
      h5_group_id = 0
      h5_prev_group_id = 0
      h5_current_time = -999.0


#endif
      END

      SUBROUTINE HDF5_CLOSE
*
*
*       Close the currently opened HDF5 handle.
*       No argument needed. This subroutine share the common block with
*       HDF5_INIT(). If the current HDF5 indicated by the h5_file_id is
*       closed, this subroutine will do nothing.
*       --------------------------
*
*
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
*       SAVE Block
      INTEGER(HID_T) :: h5_file_id, h5_step, h5_group_id
      INTEGER(HID_T), DIMENSION(1) :: h5_dset_ids(256)
      INTEGER(HID_T) :: h5_vec_len
      REAL*8 :: h5_current_time
      CHARACTER(LEN=64) :: h5_file_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER*4 :: ERROR

      IF (h5_file_id .GT. 0) THEN
        CALL h5fclose_f(h5_file_id, ERROR)
        h5_group_id = 0
        DO 5 I = 1, 256
            h5_dset_ids(I) = 0
  5    CONTINUE
        h5_step = -1
        h5_group_id = 0
        h5_prev_group_id = 0
        h5_current_time = 0
      print*, 'HDF5 CloseID:',h5_file_id,' Snapshot file: ',h5_file_name
*
      ENDIF
#endif
      END


      SUBROUTINE HDF5_write_real_vector_as_dset(group_id, dset_name,
     &     vec, vec_len, offset, dset_id, original_dset_len, finalize)
#ifdef H5OUTPUT
      USE HDF5
      IMPLICIT NONE
      INTEGER(HID_T) :: group_id ! the group ID to be written upon
      CHARACTER(LEN=*), INTENT(IN) :: dset_name ! The name of the dset
      INTEGER*4 :: vec_len ! the data array
      REAL*4  :: vec(vec_len)     ! the data array
      REAL*4  :: vec_written(vec_len)   ! the data actually written
      INTEGER*4 :: offset  ! offset of vec (by default, start from 1)
      INTEGER(HID_T) :: dset_id ! IF not 0, write to that dset; otherwise create new
      INTEGER*4 :: original_dset_len ! The original dset length, if expanding
      LOGICAL :: finalize ! Finalize the dataset, no more data can be added to it
      INTEGER(HID_T) :: I ! Loop variable
      INTEGER*4 :: error
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: memspace_id
      INTEGER(HID_T) :: crp_list ! Dataset creation property identifier
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      INTEGER(HSIZE_T), DIMENSION(2) :: data_maxdims
      INTEGER(HSIZE_T), DIMENSION(2) :: data_chunkdims
      INTEGER(HSIZE_T), DIMENSION(2) :: data_start
      INTEGER(HSIZE_T), DIMENSION(2) :: data_count

      data_dims(1) = vec_len
      data_dims(2) = 1
      data_maxdims(1) = H5S_UNLIMITED_F
      data_maxdims(2) = 1
      data_chunkdims(1) = 10
      data_chunkdims(2) = 1


      IF(dset_id .EQ. 0) THEN
*     Create new dataset
         IF(finalize .EQV. .TRUE.) THEN 
            CALL h5screate_simple_f(1, data_dims, dspace_id,error)
            CALL h5dcreate_f(group_id, dset_name, H5T_NATIVE_REAL,
     &         dspace_id, dset_id, error)
         ELSE ! Finalize = .FALSE.
*       Create simple dataspace with 1D extensible dimension
            CALL h5screate_simple_f(1, data_dims, dspace_id,
     &         error, data_maxdims)
*       Modify dataset creation properties, i.e. enable chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
*           CALL h5pset_link_creation_order_f( crp_list, 
*    &   IOR(H5P_CRT_ORDER_TRACKED_F,H5P_CRT_ORDER_INDEXED_F),error)
            CALL h5pset_chunk_f(crp_list, 1, data_chunkdims, error)
            CALL h5dcreate_f(group_id, dset_name, H5T_NATIVE_REAL,
     &         dspace_id, dset_id, error, crp_list)
         ENDIF ! finalize == .TRUE.

         IF(offset .GT. 1) THEN
*     Shift the data leftward & reduce precision
            DO 10 I = 1, vec_len
               vec_written(I) = vec(I + offset - 1)
  10        CONTINUE

cccH5_DLL herr_t H5Dwrite(hid_t dset_id, hid_t mem_type_id, hid_t mem_space_id,
ccc                         hid_t file_space_id, hid_t plist_id, const void *buf);

            CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec_written,
     &         data_dims, error)
         ELSE ! offset = 1
c$$$            DO 15 I = 1, vec_len ! only reduce precision
c$$$               vec_written(I) = REAL(vec(I))
c$$$  15        CONTINUE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec,
     &         data_dims, error)
         ENDIF ! offset > 1
         IF(finalize .EQV. .TRUE.) CALL h5dclose_f(dset_id, error)
      ELSE ! dset_id > 0
*     Expand the existing dataset
          data_start(1) = original_dset_len ! Orignal dset length (offset)
          data_start(2) = 1
          data_count(1) = vec_len ! length of the new data
          data_count(2) = 1
          data_dims(1) = vec_len + original_dset_len
          data_dims(2) = 1

          CALL h5dset_extent_f(dset_id, data_dims, error)
*       Create memspace to indicate the size of the buffer, i.e. data_count
          CALL h5screate_simple_f(1, data_count, memspace_id, error)
          CALL h5dget_space_f(dset_id, dspace_id, error)
          CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,
     &       data_start, data_count, error)

         IF(offset .GT. 1) THEN
*     Shift the data leftward & reduce precision
          DO 20 I = 1, vec_len
               vec_written(I) = vec(I + offset - 1)
  20      CONTINUE
          CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec_written,
     &        data_dims, error, memspace_id, dspace_id)
         ELSE ! offset = 1
c$$$          DO 25 I = 1, vec_len ! only reduce precision
c$$$               vec_written(I) = REAL(vec(I))
c$$$  25      CONTINUE
           CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec,
     &         data_dims, error, memspace_id, dspace_id)
         ENDIF ! offset > 1
*         CALL h5sclose_f(dspace_id, error)
         IF(finalize .EQV. .TRUE.) CALL h5dclose_f(dset_id, error)
      ENDIF ! dset_id == 0
#endif
      END


      SUBROUTINE HDF5_write_integer_vector_as_dset(group_id, dset_name,
     &     vec, vec_len, offset, dset_id, original_dset_len,finalize)
#ifdef H5OUTPUT
      USE HDF5
      INTEGER(HID_T) :: group_id ! the group ID to be written upon
      CHARACTER(LEN=*), INTENT(IN) :: dset_name ! The name of the dset
      INTEGER*4 :: vec_len ! the data array
      INTEGER*4 :: offset  ! offset of vec (by default, start from 1)
      INTEGER*4 :: vec(vec_len+offset-1)     ! the data array
      INTEGER(HID_T) :: dset_id ! IF not 0, write to that dset; otherwise create new
      INTEGER*4 :: original_dset_len ! The original dset length, before expanding
      LOGICAL :: finalize ! Finalize the dataset, no more data can be added to it
      INTEGER*4  :: vec_written(vec_len+offset-1)   !if vec starts not from 1, use this
      INTEGER(HID_T) :: I ! Loop variable
      INTEGER*4 :: error
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: memspace_id
      INTEGER(HID_T) :: crp_list ! Dataset creation property identifier
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      INTEGER(HSIZE_T), DIMENSION(2) :: data_maxdims
      INTEGER(HSIZE_T), DIMENSION(2) :: data_chunkdims
      INTEGER(HSIZE_T), DIMENSION(2) :: data_start
      INTEGER(HSIZE_T), DIMENSION(2) :: data_count

      data_dims(1) = vec_len
      data_dims(2) = 1
      data_maxdims(1) = H5S_UNLIMITED_F
      data_maxdims(2) = 1
      data_chunkdims(1) = 100
      data_chunkdims(2) = 1


      IF(dset_id .EQ. 0) THEN
*     Create new dataset
         IF(finalize .EQV. .TRUE.) THEN 
            CALL h5screate_simple_f(1, data_dims, dspace_id,error)
            CALL h5dcreate_f(group_id, dset_name, H5T_NATIVE_INTEGER,
     &         dspace_id, dset_id, error)
         ELSE 
*       Create simple dataspace with 1D extensible dimension
            CALL h5screate_simple_f(1, data_dims, dspace_id,
     &         error, data_maxdims)
*       Modify dataset creation properties, i.e. enable chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
*           CALL h5pset_link_creation_order_f( crp_list,
*    &   IOR(H5P_CRT_ORDER_TRACKED_F,H5P_CRT_ORDER_INDEXED_F),error)
            CALL h5pset_chunk_f(crp_list, 1, data_chunkdims, error)
            CALL h5dcreate_f(group_id, dset_name, H5T_NATIVE_INTEGER,
     &         dspace_id, dset_id, error, crp_list)
         ENDIF


         IF(offset .GT. 1) THEN
*     Shift the data leftward
            DO 30 I = 1, vec_len
               vec_written(I) = vec(I + offset - 1)
  30        CONTINUE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vec_written,
     &         data_dims, error)
         ELSE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vec,
     &         data_dims, error)
         ENDIF
         IF(finalize .EQV. .TRUE.) CALL h5dclose_f(dset_id, error)
      ELSE  ! dset_id > 0
*     Expand the existing dataset
          data_start(1) = original_dset_len ! Orignal dset length (offset)
          data_start(2) = 1
          data_count(1) = vec_len ! length of the new data
          data_count(2) = 1
          data_dims(1) = vec_len + original_dset_len
          data_dims(2) = 1

          CALL h5dset_extent_f(dset_id, data_dims, error)
*       Create memspace to indicate the size of the buffer, i.e. data_count
          CALL h5screate_simple_f(2, data_count, memspace_id, error)
          CALL h5dget_space_f(dset_id, dspace_id, error)
          CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,
     &       data_start, data_count, error)

         IF(offset .GT. 1) THEN
*     Shift the data leftward
          DO 40 I = 1, vec_len
               vec_written(I) = vec(I + offset - 1)
  40      CONTINUE
          CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vec_written,
     &        data_dims, error, memspace_id, dspace_id)
         ELSE ! offset = 1
           CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vec,
     &         data_dims, error, memspace_id, dspace_id)
         ENDIF ! offset > 1
         CALL h5sclose_f(dspace_id, error)
         IF(finalize .EQV. .TRUE.) CALL h5dclose_f(dset_id, error)
      ENDIF
#endif
      END

      SUBROUTINE HDF5_write_attribute_scalar_real(loc_id, att_name, val)
#ifdef H5OUTPUT
      USE HDF5
      IMPLICIT NONE
      CHARACTER(LEN=*) :: att_name
      REAL*4 :: val
      INTEGER(HID_T) :: loc_id
      INTEGER(HID_T) :: attrib_space_id
      INTEGER(HID_T) :: attrib_id
      INTEGER*4 :: error

      INTEGER(8), DIMENSION(2) :: data_dims
      data_dims(1) = 1
      data_dims(2) = 1

*       Write attributes to the group
      CALL H5Screate_f(H5S_SCALAR_F, attrib_space_id, error)
      CALL H5Acreate_f(loc_id, TRIM(att_name), H5T_NATIVE_REAL,
     &     attrib_space_id, attrib_id, error, H5P_DEFAULT_F,
     &     H5P_DEFAULT_F)
      CALL H5Awrite_f(attrib_id, H5T_NATIVE_REAL, val, data_dims,
     &      error)
      CALL H5Aclose_f(attrib_id, error)
      CALL H5Sclose_f(attrib_space_id, error)

#endif
      END

      SUBROUTINE HDF5_write_attribute_simple_real(loc_id, att_name, 
     &    vec, vec_len)
#ifdef H5OUTPUT
      USE HDF5
      IMPLICIT NONE
      CHARACTER(LEN=*) :: att_name
      INTEGER*4 :: vec_len ! the data array
      REAL*4  :: vec(vec_len)     ! the data array
      INTEGER(HID_T) :: loc_id
      INTEGER(HID_T) :: attrib_space_id
      INTEGER(HID_T) :: attrib_id
      INTEGER*4 :: error

      INTEGER(8), DIMENSION(2) :: data_dims
      data_dims(1) = vec_len
      data_dims(2) = 1

*       Write attributes to the group
      CALL H5Screate_f(H5S_SIMPLE_F, attrib_space_id, error)
      CALL H5Acreate_f(loc_id, TRIM(att_name), H5T_NATIVE_REAL,
     &     attrib_space_id, attrib_id, error, H5P_DEFAULT_F,
     &     H5P_DEFAULT_F)
      CALL H5Awrite_f(attrib_id, H5T_NATIVE_REAL, vec, data_dims,
     &      error)
      CALL H5Aclose_f(attrib_id, error)
      CALL H5Sclose_f(attrib_space_id, error)

#endif
      END


      SUBROUTINE HDF5_write_attribute_scalar_integer(loc_id, 
     &                att_name, val)
#ifdef H5OUTPUT
      USE HDF5
      IMPLICIT NONE
      CHARACTER(LEN=*) :: att_name
      INTEGER*4 :: val
      INTEGER(HID_T) :: loc_id
      INTEGER(HID_T) :: attrib_space_id
      INTEGER(HID_T) :: attrib_id
      INTEGER*4 :: error

      INTEGER(8), DIMENSION(2) :: data_dims
      data_dims(1) = 1
      data_dims(2) = 1

*       Write attributes to the group
      CALL H5Screate_f(H5S_SCALAR_F, attrib_space_id, error)
      CALL H5Acreate_f(loc_id, TRIM(att_name), H5T_NATIVE_INTEGER,
     &     attrib_space_id, attrib_id, error, H5P_DEFAULT_F,
     &     H5P_DEFAULT_F)
      CALL H5Awrite_f(attrib_id, H5T_NATIVE_INTEGER, val, data_dims,
     &      error)
      CALL H5Aclose_f(attrib_id, error)
      CALL H5Sclose_f(attrib_space_id, error)

#endif
      END


      SUBROUTINE custom_update_file(TTOT,DELTAT)
#ifdef H5OUTPUT
      USE HDF5
      REAL*8 :: TTOT, DELTAT
      REAL*8 :: h5_current_time
      INTEGER(HID_T) :: h5_file_id, h5_step, h5_group_id
      INTEGER(HID_T), DIMENSION(1) :: h5_dset_ids(256)
      INTEGER(HID_T) :: h5_vec_len
      CHARACTER(LEN=20) :: TCHAR
      CHARACTER(LEN=64) :: h5_file_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time

*     Close any previously opened file (if any)
      IF (h5_file_id .GT. 0) CALL HDF5_CLOSE
      h5_file_id = 0

      call string_left(TCHAR,TTOT,DELTAT)
      h5_file_name='snap.40_'//trim(TCHAR)//'.h5part'
#else
      COMMON/BINARYOUT/ DTOUT
      REAL*8 DELTAT,DTOUT
      DTOUT = DELTAT
#endif
      
      END

*
*     Note that N_SINGLE is now NTOT = N + NPAIRS, all members and c.m. of binaries include (R.Sp. Aug. 2021)

      SUBROUTINE output_single(TTOT,N_SINGLE,KSEV,KMODE,KRED)
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
      INCLUDE 'output_single.h'
      INTEGER(HID_T) :: h5_file_id, h5_step, h5_group_id
      INTEGER(HID_T) :: h5_subgroup_id ! ID for the Step#i groups
      INTEGER(HID_T), DIMENSION(1) :: h5_dset_ids(256)
      INTEGER(HID_T) :: h5_vec_len
      REAL*8 :: h5_current_time
      REAL*8 :: TTOT
      REAL*4 :: S_TTOT
      INTEGER*4 :: N_SINGLE,KSEV,KMODE,KRED
      CHARACTER(LEN=64) :: h5_file_name
      CHARACTER(LEN=16) :: h5_step_name
      CHARACTER(LEN=20) :: h5_step_group_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER*4 :: ERROR
      INTEGER*4 :: original_vec_len
      LOGICAL :: finalize
*
*       CALL INIT if not yet done.
      IF (h5_file_id.eq.0) then
          CALL HDF5_INIT(h5_file_name)
      END IF 
*       IF file cannot be initilized, quit the subroutine
      IF (h5_file_id .EQ. 0) RETURN

      
      IF (TTOT.GT.h5_current_time) THEN
*       Needs to create new group
          h5_step = h5_step + 1
*
          WRITE (h5_step_name, 5051) h5_step
 5051  FORMAT(I4)
*
          h5_step_group_name = 'Step#' // ADJUSTL(h5_step_name)
          CALL h5gcreate_f(h5_file_id, h5_step_group_name,
     &       h5_group_id, error) ! Create the group
*       Create subgroup for systems
*         CALL h5gcreate_f(h5_group_id, 'Systems',
*    &       h5_subgroup_id, error) ! Create the subgroup
*       Close any previously opened datasets
          DO 161 I = 1, 256 ! Close any previously opened dset
            IF(h5_dset_ids(I) .GT. 0) THEN
*               CALL h5dclose_f(h5_dset_ids(I), error)
                h5_dset_ids(I) = 0
            ENDIF
  161     CONTINUE

          S_TTOT = REAL(TTOT)
*       Write group attributes (more to follow?)
          CALL HDF5_write_attribute_scalar_real(h5_group_id, 'Time',
     &        S_TTOT) ! Write attribute to the group

          CALL HDF5_write_attribute_scalar_integer(h5_group_id,
     &    'N_STAR', N_SINGLE) ! Write attribute to the group
          
          h5_current_time = TTOT
      ELSE
*       Open existing group

cc          WRITE (h5_step_name, *), h5_step
          WRITE (h5_step_name, 5052) h5_step
5052  FORMAT(I4)

          h5_step_group_name = '/Step#' // ADJUSTL(h5_step_name)
          CALL h5gopen_f(h5_file_id, h5_step_group_name, h5_group_id,
     &        ERROR)
*       Create subgroup for singles in the opened group
*         CALL h5gcreate_f(h5_group_id, 'Singles',
*    &       h5_subgroup_id, error) ! Create the subgroup
*         CALL HDF5_write_attribute_scalar_integer(h5_subgroup_id,
*    &    'N_SINGLE', N_SINGLE) ! Write attribute to the group
      ENDIF 
      original_vec_len = 0
      finalize = .TRUE.
*
      IF (N_SINGLE.LE.0) RETURN
*
*     Write datasets
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '000 Scalars',
     &     S_AS, 100, 1, h5_dset_ids(1),
     &     original_vec_len,finalize)
*     CALL HDF5_write_attribute_simple_real(h5_group_id, 'Scalars',
*    &        S_AS,100)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '001 X1',
     &     S_X1, N_SINGLE, 1, h5_dset_ids(2),
     &     original_vec_len,finalize)
*     CALL HDF5_write_attribute_simple_real(h5_group_id,
*    &     'Pos X1',S_X1,N_SINGLE)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '002 X2',
     &     S_X2, N_SINGLE, 1, h5_dset_ids(3),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '003 X3',
     &     S_X3, N_SINGLE, 1, h5_dset_ids(4),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '004 V1',
     &     S_V1, N_SINGLE, 1, h5_dset_ids(5),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '005 V2',
     &     S_V2, N_SINGLE, 1, h5_dset_ids(6),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '006 V3',
     &     S_V3, N_SINGLE, 1, h5_dset_ids(7),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '007 A1',
     &     S_F1, N_SINGLE, 1, h5_dset_ids(8),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '008 A2',
     &     S_F2, N_SINGLE, 1, h5_dset_ids(9),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '009 A3',
     &     S_F3, N_SINGLE, 1, h5_dset_ids(10),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '010 AD1',
     &     S_FD1, N_SINGLE, 1, h5_dset_ids(11),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '011 AD2',
     &     S_FD2, N_SINGLE, 1, h5_dset_ids(12),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '012 AD3',
     &     S_FD3, N_SINGLE, 1, h5_dset_ids(13),
     &     original_vec_len,finalize)
*
      IF(KRED.EQ.0) THEN
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '013 D21',
     &     S_D21, N_SINGLE, 1, h5_dset_ids(14),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '014 D22',
     &     S_D22, N_SINGLE, 1, h5_dset_ids(15),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '015 D23',
     &     S_D23, N_SINGLE, 1, h5_dset_ids(16),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '016 D31',
     &     S_D31, N_SINGLE, 1, h5_dset_ids(17),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '017 D32',
     &     S_D32, N_SINGLE, 1, h5_dset_ids(18),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '018 D33',
     &     S_D33, N_SINGLE, 1, h5_dset_ids(19),
     &     original_vec_len,finalize)
      END IF
*
*     if(rank.eq.0)print*,' out_single d21,d31=',
*    &  (s_d21(k),k=1,n_single,100),(s_d31(k),k=1,n_single,100)
*
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '019 STEP',
     &     S_STEP, N_SINGLE, 1, h5_dset_ids(20),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '020 STEPR',
     &     S_STEPR, N_SINGLE, 1, h5_dset_ids(21),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '021 T0',
     &     S_T0, N_SINGLE, 1, h5_dset_ids(22),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '022 T0R',
     &     S_T0R, N_SINGLE, 1, h5_dset_ids(23),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '023 M',
     &     S_M, N_SINGLE, 1, h5_dset_ids(24),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '024 NB-Sph',
     &     S_RSN, N_SINGLE, 1, h5_dset_ids(25),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '025 POT',
     &     S_POT, N_SINGLE, 1, h5_dset_ids(26),
     &     original_vec_len,finalize)

      
*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '026 R*',
     &        S_RS, N_SINGLE, 1, h5_dset_ids(27),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '027 L*',
     &        S_L, N_SINGLE, 1, h5_dset_ids(28),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id,
     &        '028 Teff*', S_TE, N_SINGLE, 1, h5_dset_ids(29),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '029 RC*',
     &        S_RC, N_SINGLE, 1, h5_dset_ids(30),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '030 MC*',
     &        S_MC, N_SINGLE, 1, h5_dset_ids(31),
     &        original_vec_len,finalize)
         CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &        '031 KW', NS_KW, N_SINGLE, 1, h5_dset_ids(32),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_integer_vector_as_dset(h5_group_id,'032 Name',
     &     NS_NAM, N_SINGLE, 1, h5_dset_ids(33),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &     '033 Type', NS_LAB, N_SINGLE, 1, h5_dset_ids(34),
     &     original_vec_len,finalize)
*
      CALL h5fflush_f(h5_group_id, H5F_SCOPE_LOCAL_F, error)
      CALL h5gclose_f(h5_group_id, error) 
#else
      INCLUDE 'params.h'
      INCLUDE 'output_single.h'
      COMMON/BINARYOUT/ DTOUT
      REAL*8 DTOUT,TTOT
      INTEGER N_SINGLE,KSEV,KMODE
      CHARACTER(LEN=20) :: TCHAR
      CHARACTER(LEN=64) :: filename
      call string_left(TCHAR,TTOT,DTOUT)
      filename='single.40_'//trim(TCHAR)
      IF(KMODE.EQ.1.OR.KMODE.EQ.3) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=filename)
       WRITE(40) (S_AS(K),K=1,100)
       WRITE (40) N_SINGLE, NS_NAM(1:N_SINGLE), NS_LAB(1:N_SINGLE),
     &      S_M(1:N_SINGLE), S_RSN(1:N_SINGLE),
     &      S_X1(1:N_SINGLE), S_X2(1:N_SINGLE), S_X3(1:N_SINGLE),
     &      S_V1(1:N_SINGLE), S_V2(1:N_SINGLE), S_V3(1:N_SINGLE), 
     &      S_POT(1:N_SINGLE), S_F1(1:N_SINGLE), S_F2(1:N_SINGLE),
     &      S_POT(1:N_SINGLE), S_F1(1:N_SINGLE), S_F2(1:N_SINGLE),
     &      S_F3(1:N_SINGLE), S_FD1(1:N_SINGLE), S_FD2(1:N_SINGLE),
     &      S_FD3(1:N_SINGLE), S_STEP(1:N_SINGLE), S_STEPR(1:N_SINGLE),
     &      S_T0(1:N_SINGLE), S_T0R(1:N_SINGLE)
       IF(KRED.EQ.0) THEN
           WRITE (40)
     &     S_D21(1:N_SINGLE), S_D22(1:N_SINGLE), S_D23(1:N_SINGLE),
     &     S_D31(1:N_SINGLE), S_D32(1:N_SINGLE), S_D33(1:N_SINGLE)
       END IF
       IF(KSEV.GT.0) THEN
         WRITE (40) S_RS(1:N_SINGLE), S_L(1:N_SINGLE), S_TE(1:N_SINGLE),
     &         S_RC(1:N_SINGLE), S_MC(1:N_SINGLE), NS_KW(1:N_SINGLE)
       END IF
      ELSE IF(KMODE.EQ.2.OR.KMODE.EQ.4) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='FORMATTED',FILE=filename)
       WRITE(40,*) '## N_STAR ',N_SINGLE,' Time[NB]',TTOT
       WRITE(40,*) (S_AS(K),K=1,100)
       IF(KSEV.GT.0) THEN
          IF(KRED.EQ.1) THEN
            DO K=1,N_SINGLE
             WRITE (40,*) NS_NAM(K), NS_LAB(K), S_M(K), S_RSN(K),
     &            S_X1(K), S_X2(K), S_X3(K),
     &            S_V1(K), S_V2(K), S_V3(K), S_POT(K),
     &            S_F1(K), S_F2(K), S_F3(K), 
     &            S_FD1(K),S_FD2(K), S_FD3(K),
     &            S_STEP(K), S_STEPR(K), S_T0(K), S_T0R(K),
     &            S_RS(K), S_L(K), S_TE(K),
     &            S_MC(K), S_RC(K), NS_KW(K)
            END DO
          ELSE
            DO K=1,N_SINGLE
             WRITE (40,*) NS_NAM(K), NS_LAB(K), S_M(K), S_RSN(K),
     &            S_X1(K), S_X2(K), S_X3(K),
     &            S_V1(K), S_V2(K), S_V3(K), S_POT(K),
     &            S_F1(K), S_F2(K), S_F3(K),
     &            S_FD1(K),S_FD2(K), S_FD3(K),
     &            S_D21(K), S_D22(K), S_D23(K),
     &            S_D31(K), S_D32(K), S_D33(K),
     &            S_STEP(K), S_STEPR(K), S_T0(K), S_T0R(K),
     &            S_RS(K), S_L(K), S_TE(K),
     &            S_MC(K), S_RC(K), NS_KW(K)
            END DO
          END IF
       ELSE
          IF(KRED.EQ.1) THEN
            DO K=1,N_SINGLE
             WRITE (40,*) NS_NAM(K), NS_LAB(K), S_M(K), S_RSN(K),
     &            S_X1(K), S_X2(K), S_X3(K),
     &            S_V1(K), S_V2(K), S_V3(K), S_POT(K),
     &            S_F1(K), S_F2(K), S_F3(K),
     &            S_FD1(K),S_FD2(K), S_FD3(K),
     &            S_STEP(K), S_STEPR(K), S_T0(K), S_T0R(K)
            END DO
          ELSE
            DO K=1,N_SINGLE
             WRITE (40,*) NS_NAM(K), NS_LAB(K), S_M(K), S_RSN(K),
     &            S_X1(K), S_X2(K), S_X3(K),
     &            S_V1(K), S_V2(K), S_V3(K), S_POT(K),
     &            S_F1(K), S_F2(K), S_F3(K),
     &            S_FD1(K),S_FD2(K), S_FD3(K),
     &            S_D21(K), S_D22(K), S_D23(K),
     &            S_D31(K), S_D32(K), S_D33(K),
     &            S_STEP(K), S_STEPR(K), S_T0(K), S_T0R(K)
            END DO
          END IF
       END IF
      END IF
      CLOSE(40)
#endif
      END





      SUBROUTINE output_binary(TTOT,N_BINARY,KSEV,KMODE,KRED)
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
      INCLUDE 'output_binary.h'
      INTEGER(HID_T) :: h5_file_id, h5_step, h5_group_id
      INTEGER(HID_T), DIMENSION(1) :: h5_dset_ids(256)
      INTEGER(HID_T) :: h5_subgroup_id ! ID for the Step#i groups
      INTEGER(HID_T) :: h5_vec_len
      REAL*8 :: h5_current_time,TTOT
      REAL*4 :: S_TTOT
      INTEGER :: N_BINARY,KSEV,KMODE,KRED
      CHARACTER(LEN=64) :: h5_file_name
      CHARACTER(LEN=64) :: h5_step_name
      CHARACTER(LEN=20) :: h5_step_group_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER*4 :: ERROR
      INTEGER*4 :: original_vec_len
      LOGICAL :: finalize ! if TRUE, finalize the datasets
*
*       CALL INIT if not yet done.
      IF (h5_file_id.eq.0) then
          CALL HDF5_INIT(h5_file_name)
      END IF 
*       IF file cannot be initilized, quit the subroutine
      IF (h5_file_id .EQ. 0) RETURN

      IF (TTOT.GT.h5_current_time) THEN
*       Needs to create new group
          h5_step = h5_step + 1

cc          WRITE (h5_step_name, *), h5_step
          WRITE (h5_step_name, 5053) h5_step
5053  FORMAT(I4)

          h5_step_group_name = 'Step#' // ADJUSTL(h5_step_name)
          CALL h5gcreate_f(h5_file_id, h5_step_group_name,
     &       h5_group_id, error) ! Create the group

*       Create subgroup for binaries
*         CALL h5gcreate_f(h5_group_id, 'Binaries',
*    &       h5_subgroup_id, error) ! Create the subgroup

*       Close any previously opened datasets
*     Close any previously opened dset
          DO 361 I = 1, 256     
             IF(h5_dset_ids(I) .GT. 0) THEN
*               CALL h5dclose_f(h5_dset_ids(I), error)
                h5_dset_ids(I) = 0
             ENDIF
  361     CONTINUE
          original_vec_len = 0
          finalize = .TRUE.
          S_TTOT = REAL(TTOT)
*       Write group attributes (more to follow?)
          CALL HDF5_write_attribute_scalar_real(h5_group_id, 'Time',
     &        S_TTOT) ! Write attribute to the group

*         CALL HDF5_write_attribute_scalar_integer(h5_subgroup_id,
*    &    'N_BINARY', N_BINARY) ! Write attribute to the subgroup

          h5_current_time = TTOT
      ELSE
*       Open existing group

cc          WRITE (h5_step_name, *), h5_step
          WRITE (h5_step_name, 5054) h5_step
5054  FORMAT(I4)

          h5_step_group_name = '/Step#' // ADJUSTL(h5_step_name)
          CALL h5gopen_f(h5_file_id, h5_step_group_name, h5_group_id,
     &        ERROR)

*       Create subgroup for binaries in the opened group
*         CALL h5gcreate_f(h5_group_id, 'Binaries',
*    &       h5_subgroup_id, error) ! Create the subgroup
*         CALL HDF5_write_attribute_scalar_integer(h5_subgroup_id,
*    &    'N_BINARY', N_BINARY) ! Write attribute to the group
          original_vec_len = 0
          finalize = .TRUE.
      ENDIF

      IF (N_BINARY.LE.0) RETURN
      
*       Write datasets
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '101 Bin cm X1',
     &     B_XC1, N_BINARY, 1, h5_dset_ids(64),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '102 Bin cm X2',
     &     B_XC2, N_BINARY, 1, h5_dset_ids(65),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '103 Bin cm X3',
     &     B_XC3, N_BINARY, 1, h5_dset_ids(66),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '104 Bin cm V1',
     &     B_VC1, N_BINARY, 1, h5_dset_ids(67),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '105 Bin cm V2',
     &     B_VC2, N_BINARY, 1, h5_dset_ids(68),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '106 Bin cm V3',
     &     B_VC3, N_BINARY, 1, h5_dset_ids(69),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '107 Bin cm A1',
     &     B_FC1, N_BINARY, 1, h5_dset_ids(70),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '108 Bin cm A2',
     &     B_FC2, N_BINARY, 1, h5_dset_ids(71),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '109 Bin cm A3',
     &     B_FC3, N_BINARY, 1, h5_dset_ids(72),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '110 Bin cm AD1',
     &     B_FDC1, N_BINARY, 1, h5_dset_ids(73),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '111 Bin cm AD2',
     &     B_FDC2, N_BINARY, 1, h5_dset_ids(74),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '112 Bin cm AD3',
     &     B_FDC3, N_BINARY, 1, h5_dset_ids(75),
     &     original_vec_len,finalize)

      IF(KRED.EQ.0) THEN
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '113 Bin cm D21',
     &     B_D21, N_BINARY, 1, h5_dset_ids(76),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '114 Bin cm D22',
     &     B_D22, N_BINARY, 1, h5_dset_ids(77),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '115 Bin cm D23',
     &     B_D23, N_BINARY, 1, h5_dset_ids(78),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '116 Bin cm D31',
     &     B_D31, N_BINARY, 1, h5_dset_ids(79),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '117 Bin cm D32',
     &     B_D32, N_BINARY, 1, h5_dset_ids(80),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '118 Bin cm D33',
     &     B_D33, N_BINARY, 1, h5_dset_ids(81),
     &     original_vec_len,finalize)
      END IF

      CALL HDF5_write_real_vector_as_dset(h5_group_id, 
     &     '119 Bin cm STEP', B_STEP, N_BINARY, 1, h5_dset_ids(82),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, 
     &     '120 Bin cm STEPR', B_STEPR, N_BINARY, 1, h5_dset_ids(83),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '121 Bin cm T0',
     &     B_T0, N_BINARY, 1, h5_dset_ids(84),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '122 Bin cm T0R',
     &     B_T0R, N_BINARY, 1, h5_dset_ids(85),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '123 Bin M1*',
     &     B_M1, N_BINARY, 1, h5_dset_ids(86),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '124 Bin M2*',
     &     B_M2, N_BINARY, 1, h5_dset_ids(87),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id,'125 Bin rel X1',
     &     B_XR1, N_BINARY, 1, h5_dset_ids(88),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'126 Bin rel X2',
     &     B_XR2, N_BINARY, 1, h5_dset_ids(89),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'127 Bin rel X3',
     &     B_XR3, N_BINARY, 1, h5_dset_ids(90),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id,'128 Bin rel V1',
     &     B_VR1, N_BINARY, 1, h5_dset_ids(91),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'129 Bin rel V2',
     &     B_VR2, N_BINARY, 1, h5_dset_ids(92),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'130 Bin rel V3',
     &     B_VR3, N_BINARY, 1, h5_dset_ids(93),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id,'131 Bin rel A1',
     &     B_FR1, N_BINARY, 1, h5_dset_ids(94),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'132 Bin rel A2',
     &     B_FR2, N_BINARY, 1, h5_dset_ids(95),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'133 Bin rel A3',
     &     B_FR3, N_BINARY, 1, h5_dset_ids(96),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id,'134 Bin rel AD1',
     &     B_FDR1, N_BINARY, 1, h5_dset_ids(97),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'135 Bin rel AD2',
     &     B_FDR2, N_BINARY, 1, h5_dset_ids(98),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'136 Bin rel AD3',
     &     B_FDR3, N_BINARY, 1, h5_dset_ids(99),
     &     original_vec_len,finalize)

      IF(KRED.EQ.0) THEN
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'137 Bin rel D21',
     &     B_D21R, N_BINARY, 1, h5_dset_ids(100),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'138 Bin rel D22',
     &     B_D22R, N_BINARY, 1, h5_dset_ids(101),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'139 Bin rel D23',
     &     B_D23R, N_BINARY, 1, h5_dset_ids(102),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id,'140 Bin rel D31',
     &     B_D31R, N_BINARY, 1, h5_dset_ids(103),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'141 Bin rel D32',
     &     B_D32R, N_BINARY, 1, h5_dset_ids(104),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id,'142 Bin rel D33',
     &     B_D33R, N_BINARY, 1, h5_dset_ids(105),
     &     original_vec_len,finalize)
      END IF

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '143 Bin POT',
     &     B_POT, N_BINARY, 1, h5_dset_ids(106),
     &     original_vec_len,finalize)
      
*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_real_vector_as_dset(h5_group_id,'144 Bin RS1*',
     &        B_RS1, N_BINARY, 1, h5_dset_ids(107),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '145 Bin L1*',
     &        B_L1, N_BINARY, 1, h5_dset_ids(108),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id,
     &        '146 Bin Teff1*', B_TE1, N_BINARY, 1, h5_dset_ids(109),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id,'147 Bin RS2*',
     &        B_RS2, N_BINARY, 1, h5_dset_ids(110),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '148 Bin L2*',
     &        B_L2, N_BINARY, 1, h5_dset_ids(111),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id,
     &        '149 Bin Teff2*', B_TE2, N_BINARY, 1, h5_dset_ids(112),
     &        original_vec_len,finalize)

         CALL HDF5_write_real_vector_as_dset(h5_group_id,'150 Bin RC1*',
     &        B_RC1, N_BINARY, 1, h5_dset_ids(113),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id,'151 Bin MC1*',
     &        B_MC1, N_BINARY, 1, h5_dset_ids(114),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id,'152 Bin RC2*',
     &        B_RC2, N_BINARY, 1, h5_dset_ids(115),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id,'153 Bin MC2*',
     &        B_MC2, N_BINARY, 1, h5_dset_ids(116),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '154 Bin A[au]',
     &     B_A, N_BINARY, 1, h5_dset_ids(117),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '155 Bin ECC',
     &     B_ECC, N_BINARY, 1, h5_dset_ids(118),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '156 Bin P[d]',
     &     B_P, N_BINARY, 1, h5_dset_ids(119),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '157 Bin G',
     &     B_G, N_BINARY, 1, h5_dset_ids(120),
     &     original_vec_len,finalize)
      
*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &        '158 Bin KW1', NB_KW1, N_BINARY, 1, h5_dset_ids(121),
     &        original_vec_len,finalize)
         CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &        '159 Bin KW2', NB_KW2, N_BINARY, 1, h5_dset_ids(122),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &     '160 Bin cm KW', NB_KWC, N_BINARY, 1, h5_dset_ids(123),
     &     original_vec_len,finalize)

      CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &     '161 Bin Name1', NB_NAM1, N_BINARY, 1, h5_dset_ids(124),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &     '162 Bin Name2', NB_NAM2, N_BINARY, 1, h5_dset_ids(125),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &     '163 Bin cm Name', NB_NAMC, N_BINARY, 1, h5_dset_ids(126),
     &     original_vec_len,finalize)

      CALL h5fflush_f(h5_group_id, H5F_SCOPE_LOCAL_F, error)
      CALL h5gclose_f(h5_group_id, error) 
#else
      INCLUDE 'params.h'
      INCLUDE 'output_binary.h'
*     Binary
      COMMON/BINARYOUT/ DTOUT
      REAL*8 DTOUT,TTOT

      INTEGER N_BINARY,KSEV,KMODE
      CHARACTER(LEN=64) :: filename
      CHARACTER(LEN=20) :: TCHAR

      IF (N_BINARY.EQ.0) RETURN
      
      call string_left(TCHAR,TTOT,DTOUT)
      filename='binary.40_'//trim(TCHAR)
      IF(KMODE.EQ.1.OR.KMODE.EQ.3) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=filename)
       WRITE (40) N_BINARY, NB_NAM1(1:N_BINARY), NB_NAM2(1:N_BINARY),
     &     NB_NAMC(1:N_BINARY), B_M1(1:N_BINARY), B_M2(1:N_BINARY),
     &     B_XC1(1:N_BINARY), B_XC2(1:N_BINARY), B_XC3(1:N_BINARY),
     &     B_VC1(1:N_BINARY), B_VC2(1:N_BINARY), B_VC3(1:N_BINARY), 
     &     B_FC1(1:N_BINARY), B_FC2(1:N_BINARY), B_FC3(1:N_BINARY),
     &     B_FDC1(1:N_BINARY), B_FDC2(1:N_BINARY), B_FDC3(1:N_BINARY),
     &     B_STEP(1:N_BINARY), B_STEPR(1:N_BINARY), B_T0(1:N_BINARY),
     &     B_T0R(1:N_BINARY),
     &     B_XR1(1:N_BINARY), B_XR2(1:N_BINARY), B_XR3(1:N_BINARY), 
     &     B_VR1(1:N_BINARY), B_VR2(1:N_BINARY), B_VR3(1:N_BINARY), 
     &     B_FR1(1:N_BINARY), B_FR2(1:N_BINARY), B_FR3(1:N_BINARY),
     &     B_FDR1(1:N_BINARY), B_FDR2(1:N_BINARY), B_FDR3(1:N_BINARY),
     &     B_D21R(1:N_BINARY), B_D22R(1:N_BINARY), B_D23R(1:N_BINARY),
     &     B_D31R(1:N_BINARY), B_D32R(1:N_BINARY), B_D33R(1:N_BINARY),
     &     B_POT(1:N_BINARY),
     &     B_A(1:N_BINARY), B_ECC(1:N_BINARY), B_P(1:N_BINARY),
     &     B_G(1:N_BINARY)
       IF(KRED.EQ.0) THEN
           WRITE (40) 
     &     B_D21R(1:N_BINARY), B_D22R(1:N_BINARY), B_D23R(1:N_BINARY),
     &     B_D31R(1:N_BINARY), B_D32R(1:N_BINARY), B_D33R(1:N_BINARY)
       END IF
       IF(KSEV.GT.0) THEN
         WRITE (40) B_RS1(1:N_BINARY), B_RS2(1:N_BINARY),
     &        B_L1(1:N_BINARY), B_L2(1:N_BINARY),
     &        B_TE1(1:N_BINARY), B_TE2(1:N_BINARY),
     &        B_MC1(1:N_BINARY), B_MC2(1:N_BINARY),
     &        B_RC1(1:N_BINARY), B_RC2(1:N_BINARY), 
     &        NB_KW1(1:N_BINARY), NB_KW2(1:N_BINARY), NB_KWC(1:N_BINARY)
       END IF
      ELSE IF(KMODE.EQ.2.OR.KMODE.EQ.4) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='FORMATTED',FILE=filename)
       WRITE (40,*) '## N_BINARY ',N_BINARY, ' Time[NB]', TTOT
       IF(KSEV.GT.0) THEN
          IF(KRED.EQ.1) THEN
             DO K=1,N_BINARY
               WRITE (40,*) NB_NAM1(K), NB_NAM2(K),
     &         NB_NAMC(K), B_M1(K), B_M2(K),
     &         B_XC1(K), B_XC2(K), B_XC3(K),
     &         B_VC1(K), B_VC2(K), B_VC3(K), 
     &         B_FC1(K), B_FC2(K), B_FC3(K),
     &         B_FDC1(K), B_FDC2(K), B_FDC3(K),
     &         B_STEP(K), B_STEPR(K), B_T0(K), B_T0R(K),
     &         B_XR1(K), B_XR2(K), B_XR3(K), 
     &         B_VR1(K), B_VR2(K), B_VR3(K), B_POT(K),
     &         B_FR1(K), B_FR2(K), B_FR3(K),
     &         B_FDR1(K), B_FDR2(K), B_FDR3(K),
     &         B_A(K), B_ECC(K), B_P(K), B_G(K),
     &         B_RS1(K), B_RS2(K),
     &         B_L1(K), B_L2(K),
     &         B_TE1(K), B_TE2(K),
     &         B_MC1(K), B_MC2(K),
     &         B_RC1(K), B_RC2(K), 
     &         NB_KW1(K), NB_KW2(K), NB_KWC(K)
             END DO
          ELSE 
             DO K=1,N_BINARY
             WRITE (40,*) NB_NAM1(K), NB_NAM2(K),
     &         NB_NAMC(K), B_M1(K), B_M2(K),
     &         B_XC1(K), B_XC2(K), B_XC3(K),
     &         B_VC1(K), B_VC2(K), B_VC3(K),
     &         B_FC1(K), B_FC2(K), B_FC3(K),
     &         B_FDC1(K), B_FDC2(K), B_FDC3(K),
     &         B_STEP(K), B_STEPR(K), B_T0(K), B_T0R(K),
     &         B_XR1(K), B_XR2(K), B_XR3(K),
     &         B_VR1(K), B_VR2(K), B_VR3(K), B_POT(K),
     &         B_FR1(K), B_FR2(K), B_FR3(K),
     &         B_FDR1(K), B_FDR2(K), B_FDR3(K),
     &         B_D21R(K), B_D22R(K), B_D23R(K),
     &         B_D31R(K), B_D32R(K), B_D33R(K),
     &         B_A(K), B_ECC(K), B_P(K), B_G(K),
     &         B_RS1(K), B_RS2(K),
     &         B_L1(K), B_L2(K),
     &         B_TE1(K), B_TE2(K),
     &         B_MC1(K), B_MC2(K),
     &         B_RC1(K), B_RC2(K),
     &         NB_KW1(K), NB_KW2(K), NB_KWC(K)
             END DO
          END IF
       ELSE
          IF(KRED.EQ.1) THEN
              DO K=1,N_BINARY
                WRITE (40,*) NB_NAM1(K), NB_NAM2(K),
     &          NB_NAMC(K), B_M1(K), B_M2(K),
     &          B_XC1(K), B_XC2(K), B_XC3(K),
     &          B_VC1(K), B_VC2(K), B_VC3(K), 
     &          B_FC1(K), B_FC2(K), B_FC3(K),
     &          B_FDC1(K), B_FDC2(K), B_FDC3(K),
     &          B_STEP(K), B_STEPR(K), B_T0(K), B_T0R(K),
     &          B_XR1(K), B_XR2(K), B_XR3(K), 
     &          B_VR1(K), B_VR2(K), B_VR3(K), B_POT(K),
     &          B_FR1(K), B_FR2(K), B_FR3(K),
     &          B_FDR1(K), B_FDR2(K), B_FDR3(K),
     &          B_A(K), B_ECC(K), B_P(K)
              END DO
          ELSE
              DO K=1,N_BINARY
                WRITE (40,*) NB_NAM1(K), NB_NAM2(K),
     &          NB_NAMC(K), B_M1(K), B_M2(K),
     &          B_XC1(K), B_XC2(K), B_XC3(K),
     &          B_VC1(K), B_VC2(K), B_VC3(K),
     &          B_FC1(K), B_FC2(K), B_FC3(K),
     &          B_FDC1(K), B_FDC2(K), B_FDC3(K),
     &          B_STEP(K), B_STEPR(K), B_T0(K), B_T0R(K),
     &          B_XR1(K), B_XR2(K), B_XR3(K),
     &          B_VR1(K), B_VR2(K), B_VR3(K), B_POT(K),
     &          B_FR1(K), B_FR2(K), B_FR3(K),
     &          B_FDR1(K), B_FDR2(K), B_FDR3(K),
     &          B_D21R(K), B_D22R(K), B_D23R(K),
     &          B_D31R(K), B_D32R(K), B_D33R(K),
     &          B_A(K), B_ECC(K), B_P(K)
              END DO
         END IF
        END IF
       END IF

         
      CLOSE(40)
      
#endif
      END





      SUBROUTINE output_merger(TTOT,N_MERGER,KSEV,KMODE,KRED)
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
      INCLUDE 'output_merger.h'
      INTEGER(HID_T) :: h5_file_id, h5_step, h5_group_id
      INTEGER(HID_T), DIMENSION(1) :: h5_dset_ids(256)
      INTEGER(HID_T) :: h5_subgroup_id
      INTEGER(HID_T) :: h5_vec_len
      REAL*8 :: h5_current_time,TTOT
      REAL*4 :: S_TTOT
      INTEGER :: N_MERGER,KSEV,KMODE,KRED
      CHARACTER(LEN=64) :: h5_file_name
      CHARACTER(LEN=64) :: h5_step_name
      CHARACTER(LEN=20) :: h5_step_group_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER*4 :: ERROR
      INTEGER*4 :: original_vec_len
      LOGICAL :: finalize

C*     Merger
*       CALL INIT if not yet done.
      IF (h5_file_id.eq.0) then
          CALL HDF5_INIT(h5_file_name)
      END IF 
*       IF file cannot be initilized, quit the subroutine
      IF (h5_file_id .EQ. 0) RETURN

      
      IF (TTOT.GT.h5_current_time) THEN
*       Needs to create new group
           h5_step = h5_step + 1 

cc          WRITE (h5_step_name, *), h5_step
          WRITE (h5_step_name, 5055) h5_step
5055  FORMAT(I4)

          h5_step_group_name = 'Step#' // ADJUSTL(h5_step_name)
          CALL h5gcreate_f(h5_file_id, h5_step_group_name,
     &       h5_group_id, error) ! Create the group

*       Create subgroup for mergers
*         CALL h5gcreate_f(h5_group_id, 'Mergers',
*    &       h5_subgroup_id, error) ! Create the subgroup

*       Close any previously opened datasets
          DO 361 I = 1, 256     ! Close any previously opened dset
             IF(h5_dset_ids(I) .GT. 0) THEN
*     CALL h5dclose_f(h5_dset_ids(I), error)
                h5_dset_ids(I) = 0
             ENDIF
  361    CONTINUE
         original_vec_len = 0
         finalize = .TRUE.
         S_TTOT = REAL(TTOT)
*       Write group attributes (more to follow?)
          CALL HDF5_write_attribute_scalar_real(h5_group_id, 'Time',
     &        S_TTOT) ! Write attribute to the group

*         CALL HDF5_write_attribute_scalar_integer(h5_subgroup_id,
*    &    'N_MERGER', N_MERGER) ! Write attribute to the group

          h5_current_time = TTOT

      ELSE
*       Open existing group

cc          WRITE (h5_step_name, *), h5_step
          WRITE (h5_step_name, 5056) h5_step
5056  FORMAT(I4)


          h5_step_group_name = '/Step#' // ADJUSTL(h5_step_name)
          CALL h5gopen_f(h5_file_id, h5_step_group_name, h5_group_id,
     &        ERROR)

*       Create subgroup for mergers in the opened group
*         CALL h5gcreate_f(h5_group_id, 'Mergers',
*    &       h5_subgroup_id, error)  ! Create the subgroup
*         CALL HDF5_write_attribute_scalar_integer(h5_subgroup_id,
*    &    'N_MERGER', N_MERGER) ! Write attribute to the group
          original_vec_len = 0
          finalize = .TRUE.
      ENDIF 

      IF (N_MERGER.LE.0) RETURN

*       Write datasets
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '201 Mer XC1',
     &     M_XC1, N_MERGER, 1, h5_dset_ids(182),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '202 Mer XC2',
     &     M_XC2, N_MERGER, 1, h5_dset_ids(183),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '203 Mer XC3',
     &     M_XC3, N_MERGER, 1, h5_dset_ids(184),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '204 Mer VC1',
     &     M_VC1, N_MERGER, 1, h5_dset_ids(185),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '205 Mer VC2',
     &     M_VC2, N_MERGER, 1, h5_dset_ids(186),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '206 Mer VC3',
     &     M_VC3, N_MERGER, 1, h5_dset_ids(187),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '207 Mer M1',
     &     M_M1, N_MERGER, 1, h5_dset_ids(188),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '208 Mer M2',
     &     M_M2, N_MERGER, 1, h5_dset_ids(189),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '209 Mer M3',
     &     M_M3, N_MERGER, 1, h5_dset_ids(190),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '210 Mer XR01',
     &     M_XR01, N_MERGER, 1, h5_dset_ids(191),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '211 Mer XR02',
     &     M_XR02, N_MERGER, 1, h5_dset_ids(192),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '212 Mer XR03',
     &     M_XR03, N_MERGER, 1, h5_dset_ids(193),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '213 Mer VR01',
     &     M_VR01, N_MERGER, 1, h5_dset_ids(194),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '214 Mer VR02',
     &     M_VR02, N_MERGER, 1, h5_dset_ids(195),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '215 Mer VR03',
     &     M_VR03, N_MERGER, 1, h5_dset_ids(196),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '216 Mer XR11',
     &     M_XR11, N_MERGER, 1, h5_dset_ids(197),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '217 Mer XR12',
     &     M_XR12, N_MERGER, 1, h5_dset_ids(198),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '218 Mer XR13',
     &     M_XR13, N_MERGER, 1, h5_dset_ids(199),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '219 Mer VR11',
     &     M_VR11, N_MERGER, 1, h5_dset_ids(200),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '220 Mer VR12',
     &     M_VR12, N_MERGER, 1, h5_dset_ids(201),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '221 Mer VR13',
     &     M_VR13, N_MERGER, 1, h5_dset_ids(202),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '222 Mer POT',
     &     M_POT, N_MERGER, 1, h5_dset_ids(203),
     &     original_vec_len,finalize)

*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '223 Mer RS1',
     &        M_RS1, N_MERGER, 1, h5_dset_ids(204),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '224 Mer L1',
     &        M_L1, N_MERGER, 1, h5_dset_ids(205),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '225 Mer TE1',
     &        M_TE1, N_MERGER, 1, h5_dset_ids(206),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '226 Mer RS2',
     &        M_RS2, N_MERGER, 1, h5_dset_ids(207),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '227 Mer L2',
     &        M_L2, N_MERGER, 1, h5_dset_ids(208),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '228 Mer TE2',
     &        M_TE2, N_MERGER, 1, h5_dset_ids(209),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '229 Mer RS3',
     &        M_RS3, N_MERGER, 1, h5_dset_ids(210),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '230 Mer L3',
     &        M_L3, N_MERGER, 1, h5_dset_ids(211),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '231 Mer TE3',
     &        M_TE3, N_MERGER, 1, h5_dset_ids(212),
     &        original_vec_len,finalize)

         CALL HDF5_write_real_vector_as_dset(h5_group_id, '232 Mer RC1',
     &        M_RC1, N_MERGER, 1, h5_dset_ids(213),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '233 Mer MC1',
     &        M_MC1, N_MERGER, 1, h5_dset_ids(214),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '234 Mer RC2',
     &        M_RC2, N_MERGER, 1, h5_dset_ids(215),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '235 Mer MC2',
     &        M_MC2, N_MERGER, 1, h5_dset_ids(216),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '236 Mer RC3',
     &        M_RC3, N_MERGER, 1, h5_dset_ids(217),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, '237 Mer MC3',
     &        M_MC3, N_MERGER, 1, h5_dset_ids(218),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_real_vector_as_dset(h5_group_id, '238 Mer A0[au]',
     &     M_A0, N_MERGER, 1, h5_dset_ids(219),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '239 Mer ECC0',
     &     M_ECC0, N_MERGER, 1, h5_dset_ids(220),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '240 Mer P0[d]',
     &     M_P0, N_MERGER, 1, h5_dset_ids(221),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '241 Mer A1[au]',
     &     M_A1, N_MERGER, 1, h5_dset_ids(222),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '242 Mer ECC1',
     &     M_ECC1, N_MERGER, 1, h5_dset_ids(223),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, '243 Mer P1[d]',
     &     M_P1, N_MERGER, 1, h5_dset_ids(224),
     &     original_vec_len,finalize)

*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &        '244 Mer KW1', NM_KW1, N_MERGER, 1, h5_dset_ids(225),
     &        original_vec_len,finalize)
         CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &        '245 Mer KW2', NM_KW2, N_MERGER, 1, h5_dset_ids(226),
     &        original_vec_len,finalize)
         CALL HDF5_write_integer_vector_as_dset(h5_group_id,
     &        '246 Mer KW3', NM_KW3, N_MERGER, 1, h5_dset_ids(227),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_integer_vector_as_dset(h5_group_id, '247 Mer KWC',
     &     NM_KWC, N_MERGER, 1, h5_dset_ids(228),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_group_id,'248 Mer NAM1',
     &     NM_NAM1, N_MERGER, 1, h5_dset_ids(229),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_group_id,'249 Mer NAM2',
     &     NM_NAM2, N_MERGER, 1, h5_dset_ids(230),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_group_id,'250 Mer NAM3',
     &     NM_NAM3, N_MERGER, 1, h5_dset_ids(231),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_group_id,'251 Mer NAMC',
     &     NM_NAMC, N_MERGER, 1, h5_dset_ids(232),
     &     original_vec_len,finalize)

      CALL h5fflush_f(h5_group_id, H5F_SCOPE_LOCAL_F, error)
      CALL h5gclose_f(h5_group_id, error) 

#else
      INCLUDE 'params.h'
      INCLUDE 'output_merger.h'
      COMMON/BINARYOUT/ DTOUT
      REAL*8 DTOUT,TTOT

      INTEGER N_MERGER,KSEV,KMODE
      CHARACTER(LEN=20) :: TCHAR
      CHARACTER(LEN=64) :: filename

      IF (N_MERGER.EQ.0) RETURN
      
      call string_left(TCHAR,TTOT,DTOUT)
      filename='merger.40_'//trim(TCHAR)
      IF(KMODE.EQ.1.OR.KMODE.EQ.3) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=filename)
       WRITE (40) NM_NAM1(1:N_MERGER), NM_NAM2(1:N_MERGER),
     &     NM_NAM3(1:N_MERGER), NM_NAMC(1:N_MERGER),
     &     M_M1(1:N_MERGER), M_M2(1:N_MERGER), M_M3(1:N_MERGER),
     &     M_XC1(1:N_MERGER), M_XC2(1:N_MERGER), M_XC3(1:N_MERGER),
     &     M_VC1(1:N_MERGER), M_VC2(1:N_MERGER), M_VC3(1:N_MERGER),
     &     M_XR01(1:N_MERGER), M_XR02(1:N_MERGER), M_XR03(1:N_MERGER), 
     &     M_VR01(1:N_MERGER), M_VR02(1:N_MERGER), M_VR03(1:N_MERGER), 
     &     M_XR11(1:N_MERGER), M_XR12(1:N_MERGER), M_XR13(1:N_MERGER), 
     &     M_VR11(1:N_MERGER), M_VR12(1:N_MERGER), M_VR13(1:N_MERGER),
     &     M_POT(1:N_MERGER)
       IF (KSEV.GT.0) THEN
      WRITE(40) M_RS1(1:N_MERGER), M_RS2(1:N_MERGER), M_RS3(1:N_MERGER),
     &        M_L1(1:N_MERGER), M_L2(1:N_MERGER), M_L3(1:N_MERGER), 
     &        M_TE1(1:N_MERGER), M_TE2(1:N_MERGER), M_TE3(1:N_MERGER), 
     &        M_RC1(1:N_MERGER), M_RC2(1:N_MERGER), M_RC3(1:N_MERGER),
     &        M_MC1(1:N_MERGER), M_MC2(1:N_MERGER), M_MC3(1:N_MERGER),
     &        M_A0(1:N_MERGER), M_ECC0(1:N_MERGER), M_P0(1:N_MERGER),
     &        M_A1(1:N_MERGER), M_ECC1(1:N_MERGER), M_P1(1:N_MERGER),
     &        NM_KW1(1:N_MERGER), NM_KW2(1:N_MERGER),
     &        NM_KW3(1:N_MERGER), NM_KWC(1:N_MERGER)
       END IF
      ELSE IF(KMODE.EQ.2.OR.KMODE.EQ.4) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='FORMATTED',FILE=filename)
       WRITE (40,*) '## N_MERGER ',N_MERGER, ' Time[NB]', TTOT
       IF(KSEV.GT.0) THEN
          DO K=1,N_MERGER
             WRITE (40,*) NM_NAM1(K), NM_NAM2(K),
     &            NM_NAM3(K), NM_NAMC(K),
     &            M_M1(K), M_M2(K), M_M3(K),
     &            M_XC1(K), M_XC2(K), M_XC3(K),
     &            M_VC1(K), M_VC2(K), M_VC3(K),
     &            M_XR01(K), M_XR02(K), M_XR03(K), 
     &            M_VR01(K), M_VR02(K), M_VR03(K), 
     &            M_XR11(K), M_XR12(K), M_XR13(K), 
     &            M_VR11(K), M_VR12(K), M_VR13(K), M_POT(K),
     &            M_A0(K), M_ECC0(K), M_P0(K),
     &            M_A1(K), M_ECC1(K), M_P1(K),
     &            M_RS1(K), M_RS2(K), M_RS3(K),
     &            M_L1(K), M_L2(K), M_L3(K), 
     &            M_TE1(K), M_TE2(K), M_TE3(K), 
     &            M_MC1(K), M_MC2(K), M_MC3(K),
     &            M_RC1(K), M_RC2(K), M_RC3(K),
     &            NM_KW1(K), NM_KW2(K),
     &            NM_KW3(K), NM_KWC(K)
          END DO
       ELSE
          DO K=1,N_MERGER
             WRITE (40,*) NM_NAM1(K), NM_NAM2(K),
     &            NM_NAM3(K), NM_NAMC(K),
     &            M_M1(K), M_M2(K), M_M3(K),
     &            M_XC1(K), M_XC2(K), M_XC3(K),
     &            M_VC1(K), M_VC2(K), M_VC3(K),
     &            M_XR01(K), M_XR02(K), M_XR03(K), 
     &            M_VR01(K), M_VR02(K), M_VR03(K), 
     &            M_XR11(K), M_XR12(K), M_XR13(K), 
     &            M_VR11(K), M_VR12(K), M_VR13(K), M_POT(K),
     &            M_A0(K), M_ECC0(K), M_P0(K),
     &            M_A1(K), M_ECC1(K), M_P1(K)
          END DO
       END IF
      END IF
      CLOSE(40)
      
#endif
      END
