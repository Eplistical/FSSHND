#ifndef _MATRIXOP_CONFIG_HPP
#define _MATRIXOP_CONFIG_HPP

#include <complex>

// types
#ifndef MATRIXOP_TYPEDEF

#define MATRIXOP_CONST const
#define MATRIXOP_CHAR  char
#define MATRIXOP_VOID  void
#define MATRIXOP_BOOL  bool
#define MATRIXOP_ITYPE int
#define MATRIXOP_STYPE float
#define MATRIXOP_DTYPE double
#define MATRIXOP_CTYPE std::complex<float>
#define MATRIXOP_ZTYPE std::complex<double>

#endif 

// static
#ifdef _OPENMP
    #define MATRIXOP_STATIC 
#else 
    #define MATRIXOP_STATIC static
#endif

#define MATRIXOP_STRICT_STATIC static

// Interface 
#ifndef MATRIXOP_INTERFACE
#define MATRIXOP_INTERFACE_NETLIB
#endif

// Fortran Naming Convention
#ifndef MATRIXOP_FORTRAN_CONVENTION
#define MATRIXOP_FORTRAN_UNDERSCORE
#endif 


// Interface Macro Def
#ifdef MATRIXOP_FORTRAN_UNDERSCORE


// -- BLAS LEVEL 1 -- //


// single
#define SASUM       sasum_
#define SAXPY       saxpy_
#define SCABS1      scabs1_
#define SCASUM      scasum_
#define SCNRM2      scnrm2_
#define SCOPY       scopy_
#define SDOT        sdot_
#define SDSDOT      sdsdot_
#define SNRM2       snrm2_
#define SROT        srot_
#define SROTG       srotg_
#define SROTM       srotm_
#define SROTMG      srotmg_
#define SSCAL       sscal_
#define SSWAP       sswap_
// double
#define DASUM       dasum_
#define DAXPY       daxpy_
#define DCABS1      dcabs1_
#define DCOPY       dcopy_
#define DDOT        ddot_
#define DNRM2       dnrm2_
#define DROT        drot_
#define DROTG       drotg_
#define DROTM       drotm_
#define DROTMG      drotmg_
#define DSCAL       dscal_
#define DSDOT       dsdot_
#define DSWAP       dswap_
#define DZASUM      dzasum_
#define DZNRM2      dznrm2_
// complex
#define CAXPY       caxpy_
#define CCOPY       ccopy_
#define CDOTC       cdotc_
#define CDOTU       cdotu_
#define CROTG       crotg_
#define CSCAL       cscal_
#define CSROT       csrot_
#define CSSCAL      csscal_
#define CSWAP       cswap_
// complex*16
#define ZAXPY       zaxpy_
#define ZCOPY       zcopy_
#define ZDOTC       zdotc_
#define ZDOTU       zdotu_
#define ZDROT       zdrot_
#define ZDSCAL      zdscal_
#define ZROTG       zrotg_
#define ZSCAL       zscal_
#define ZSWAP       zswap_


// -- BLAS LEVEL 1 END -- //


// -- BLAS LEVEL 2 -- //


// single
#define SGBMV		sgbmv_
#define SGEMV		sgemv_
#define SGER		sger_
#define SSBMV		ssbmv_
#define SSPMV		sspmv_
#define SSPR		sspr_
#define SSPR2		sspr2_
#define SSYMV		ssymv_
#define SSYR		ssyr_
#define SSYR2		ssyr2_
#define STBMV		stbmv_
#define STBSV		stbsv_
#define STPMV		stpmv_
#define STPSV		stpsv_
#define STRMV		strmv_
#define STRSV		strsv_
// double
#define DGBMV		dgbmv_
#define DGEMV		dgemv_
#define DGER		dger_
#define DSBMV		dsbmv_
#define DSPMV		dspmv_
#define DSPR		dspr_
#define DSPR2		dspr2_
#define DSYMV		dsymv_
#define DSYR		dsyr_
#define DSYR2		dsyr2_
#define DTBMV		dtbmv_
#define DTBSV		dtbsv_
#define DTPMV		dtpmv_
#define DTPSV		dtpsv_
#define DTRMV		dtrmv_
#define DTRSV		dtrsv_
// complex
#define CGBMV		cgbmv_
#define CGEMV		cgemv_
#define CGERC		cgerc_
#define CGERU		cgeru_
#define CHBMV		chbmv_
#define CHEMV		chemv_
#define CHER		cher_
#define CHER2		cher2_
#define CHPMV		chpmv_
#define CHPR		chpr_
#define CHPR2		chpr2_
#define CTBMV		ctbmv_
#define CTBSV		ctbsv_
#define CTPMV		ctpmv_
#define CTPSV		ctpsv_
#define CTRMV		ctrmv_
#define CTRSV		ctrsv_
// complex*16
#define ZGBMV		zgbmv_
#define ZGEMV		zgemv_
#define ZGERC		zgerc_
#define ZGERU		zgeru_
#define ZHBMV		zhbmv_
#define ZHEMV		zhemv_
#define ZHER		zher_
#define ZHER2		zher2_
#define ZHPMV		zhpmv_
#define ZHPR		zhpr_
#define ZHPR2		zhpr2_
#define ZTBMV		ztbmv_
#define ZTBSV		ztbsv_
#define ZTPMV		ztpmv_
#define ZTPSV		ztpsv_
#define ZTRMV		ztrmv_
#define ZTRSV		ztrsv_


// -- BLAS LEVEL 2 END -- //


// -- BLAS LEVEL 3 -- //


// single
#define SGEMM		sgemm_
#define SSYMM		ssymm_
#define SSYR2K		ssyr2k_
#define SSYRK		ssyrk_
#define STRMM		strmm_
#define STRSM		strsm_
// double
#define DGEMM		dgemm_
#define DSYMM		dsymm_
#define DSYR2K		dsyr2k_
#define DSYRK		dsyrk_
#define DTRMM		dtrmm_
#define DTRSM		dtrsm_
// complex
#define CGEMM		cgemm_
#define CHEMM		chemm_
#define CHER2K		cher2k_
#define CHERK		cherk_
#define CSYMM		csymm_
#define CSYR2K		csyr2k_
#define CSYRK		csyrk_
#define CTRMM		ctrmm_
#define CTRSM		ctrsm_
// complex*16
#define ZGEMM		zgemm_
#define ZHEMM		zhemm_
#define ZHER2K		zher2k_
#define ZHERK		zherk_
#define ZSYMM		zsymm_
#define ZSYR2K		zsyr2k_
#define ZSYRK		zsyrk_
#define ZTRMM		ztrmm_
#define ZTRSM		ztrsm_


// -- BLAS LEVEL 3 END -- //


// -- LAPACK -- //


// single
#define SGEES       sgees_
#define SGETRF      sgetrf_
#define SGETRI      sgetri_
#define SSYEVD      ssyevd_
#define SGEEV       sgeev_
#define SGESV       sgesv_

// double
#define DGEES       dgees_
#define DGETRF      dgetrf_
#define DGETRI      dgetri_
#define DSYEVD      dsyevd_
#define DGEEV       dgeev_
#define DGESV       dgesv_

// complex
#define CGEES       cgees_
#define CGETRF      cgetrf_
#define CGETRI      cgetri_
#define CHEEVD      cheevd_
#define CGEEV       cgeev_
#define CGESV       cgesv_

// complex*16
#define ZGEES       zgees_
#define ZGETRF      zgetrf_
#define ZGETRI      zgetri_
#define ZHEEVD      zheevd_
#define ZGEEV       zgeev_
#define ZGESV       zgesv_


// -- LAPACK END -- //


#endif // MATRIXOP_FORTRAN_UNDERSCORE


namespace matrixop {

    // types
    using VOID  = MATRIXOP_VOID;

    using CHAR  = MATRIXOP_CHAR;
    using BOOL  = MATRIXOP_BOOL;
    using ITYPE = MATRIXOP_ITYPE;
    using STYPE = MATRIXOP_STYPE;
    using DTYPE = MATRIXOP_DTYPE;
    using CTYPE = MATRIXOP_CTYPE;
    using ZTYPE = MATRIXOP_ZTYPE;

    using CNST_CHAR  = MATRIXOP_CONST CHAR;
    using CNST_BOOL  = MATRIXOP_CONST BOOL;
    using CNST_ITYPE = MATRIXOP_CONST ITYPE;
    using CNST_STYPE = MATRIXOP_CONST STYPE;
    using CNST_DTYPE = MATRIXOP_CONST DTYPE;
    using CNST_CTYPE = MATRIXOP_CONST CTYPE;
    using CNST_ZTYPE = MATRIXOP_CONST ZTYPE;


    // function types

    using SELE1_STYPE = BOOL (*)(CNST_STYPE*);
    using SELE1_DTYPE = BOOL (*)(CNST_DTYPE*);
    using SELE1_CTYPE = BOOL (*)(CNST_CTYPE*);
    using SELE1_ZTYPE = BOOL (*)(CNST_ZTYPE*);

    using SELE2_STYPE = BOOL (*)(CNST_STYPE*, CNST_STYPE*);
    using SELE2_DTYPE = BOOL (*)(CNST_DTYPE*, CNST_DTYPE*);
    using SELE2_CTYPE = BOOL (*)(CNST_CTYPE*, CNST_CTYPE*);
    using SELE2_ZTYPE = BOOL (*)(CNST_ZTYPE*, CNST_ZTYPE*);

    // BLAS/LAPACK constants
	static CNST_CHAR CHARC('C');
	static CNST_CHAR CHARE('E');
	static CNST_CHAR CHARF('F');
	static CNST_CHAR CHARL('L');
	static CNST_CHAR CHARN('N');
	static CNST_CHAR CHARS('S');
	static CNST_CHAR CHART('T');
	static CNST_CHAR CHARV('V');
	static CNST_CHAR CHARU('U');

	static CNST_ITYPE ONEI(1);
	static CNST_STYPE ONES(1.0);
	static CNST_DTYPE ONED(1.0);
	static CNST_CTYPE ONEC(1.0, 0.0);
	static CNST_ZTYPE ONEZ(1.0, 0.0);

	static CNST_ITYPE ZEROI(0);
	static CNST_STYPE ZEROS(0.0);
	static CNST_DTYPE ZEROD(0.0);
	static CNST_CTYPE ZEROC(0.0, 0.0);
	static CNST_ZTYPE ZEROZ(0.0, 0.0);

	static CNST_CTYPE IMAGIC(0.0, 1.0);
	static CNST_ZTYPE IMAGIZ(0.0, 1.0);

    // BLAS/LAPACK interfaces
	extern "C" {


        // -- BLAS LEVEL 1 -- //


        // single
        extern STYPE    SASUM		(CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*);
        extern VOID     SAXPY		(CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        extern STYPE    SCABS1		(CNST_CTYPE*);
        extern STYPE    SCASUM		(CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*);
        extern STYPE    SCNRM2		(CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*);
        extern VOID     SCOPY		(CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        extern STYPE    SDOT		(CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*);
        extern STYPE    SDSDOT		(CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*);
        extern STYPE    SNRM2		(CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*);
        extern VOID     SROT		(CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*);
        extern VOID     SROTG		(     STYPE*,      STYPE*,      STYPE*,      STYPE*);
        extern VOID     SROTM		(CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*, CNST_STYPE*);
        extern VOID     SROTMG		(     STYPE*,      STYPE*,      STYPE*, CNST_STYPE*,      STYPE*);
        extern VOID     SSCAL		(CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID     SSWAP		(CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        // double
        extern DTYPE    DASUM		(CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*);
        extern VOID     DAXPY		(CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern DTYPE    DCABS1		(CNST_ZTYPE*);
        extern VOID     DCOPY		(CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern DTYPE    DDOT		(CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*);
        extern DTYPE    DNRM2		(CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*);
        extern VOID     DROT		(CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*);
        extern VOID     DROTG		(     DTYPE*,      DTYPE*,      DTYPE*,      DTYPE*);
        extern VOID     DROTM		(CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*, CNST_DTYPE*);
        extern VOID     DROTMG		(     DTYPE*,      DTYPE*,      DTYPE*, CNST_DTYPE*,      DTYPE*);
        extern VOID     DSCAL		(CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern DTYPE    DSDOT		(CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*);
        extern VOID     DSWAP		(CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern DTYPE    DZASUM		(CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*);
        extern DTYPE    DZNRM2		(CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*);
        // complex
        extern VOID     CAXPY		(CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID     CCOPY		(CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern CTYPE    CDOTC		(CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*);
        extern CTYPE    CDOTU		(CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*);
        extern VOID     CROTG		(     CTYPE*, CNST_CTYPE*,      STYPE*,      CTYPE*);
        extern VOID     CSCAL		(CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID     CSROT		(CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*);
        extern VOID     CSSCAL		(CNST_ITYPE*, CNST_STYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID     CSWAP		(CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        // complex*16
        extern VOID     ZAXPY		(CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID     ZCOPY		(CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern ZTYPE    ZDOTC		(CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*);
        extern ZTYPE    ZDOTU		(CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*);
        extern VOID     ZDROT		(CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*);
        extern VOID     ZDSCAL		(CNST_ITYPE*, CNST_DTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID     ZROTG		(     ZTYPE*, CNST_ZTYPE*,      DTYPE*,      ZTYPE*);
        extern VOID     ZSCAL		(CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID     ZSWAP		(CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);


        // -- BLAS LEVEL 1 END -- //


        // -- BLAS LEVEL 2 -- //
        

        // single
        extern VOID		SGBMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SGEMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SGER		(CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SSBMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SSPMV		( CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SSPR		( CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*);
        extern VOID		SSPR2		( CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*);
        extern VOID		SSYMV		( CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SSYR		( CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SSYR2		( CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		STBMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		STBSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		STPMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		STPSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		STRMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		STRSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        // double
        extern VOID		DGBMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DGEMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DGER		(CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DSBMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DSPMV		( CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DSPR		( CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*);
        extern VOID		DSPR2		( CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*);
        extern VOID		DSYMV		( CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DSYR		( CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DSYR2		( CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DTBMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DTBSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DTPMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DTPSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DTRMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DTRSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        // complex
        extern VOID		CGBMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CGEMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CGERC		(CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CGERU		(CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CHBMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CHEMV		( CNST_CHAR*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CHER		( CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CHER2		( CNST_CHAR*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CHPMV		( CNST_CHAR*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CHPR		( CNST_CHAR*, CNST_ITYPE*, CNST_STYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*);
        extern VOID		CHPR2		( CNST_CHAR*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*);
        extern VOID		CTBMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CTBSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CTPMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CTPSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CTRMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CTRSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        // complex*16
        extern VOID		ZGBMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZGEMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZGERC		(CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZGERU		(CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZHBMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZHEMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZHER		( CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZHER2		( CNST_CHAR*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZHPMV		( CNST_CHAR*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZHPR		( CNST_CHAR*, CNST_ITYPE*, CNST_DTYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*);
        extern VOID		ZHPR2		( CNST_CHAR*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*);
        extern VOID		ZTBMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZTBSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZTPMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZTPSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*,      ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZTRMV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZTRSV		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);


        // -- BLAS LEVEL 2 END -- //


        // -- BLAS LEVEL 3 -- //


        // single
        extern VOID		SGEMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SSYMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SSYR2K		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		SSYRK		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*, CNST_STYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		STRMM		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        extern VOID		STRSM		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*);
        // double
        extern VOID		DGEMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DSYMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DSYR2K		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DSYRK		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*, CNST_DTYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DTRMM		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        extern VOID		DTRSM		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*);
        // complex
        extern VOID		CGEMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CHEMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CHER2K		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_STYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CHERK		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_STYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_STYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CSYMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CSYR2K		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CSYRK		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*, CNST_CTYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CTRMM		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        extern VOID		CTRSM		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_CTYPE*, CNST_CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*);
        // complex*16
        extern VOID		ZGEMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZHEMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZHER2K		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_DTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZHERK		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_DTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_DTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZSYMM		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZSYR2K		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZSYRK		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*, CNST_ZTYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZTRMM		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);
        extern VOID		ZTRSM		( CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*, CNST_ITYPE*, CNST_ZTYPE*, CNST_ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*);


        // -- BLAS LEVEL 3 END -- //



        // -- lapack -- //

        // gees
        extern VOID		SGEES		( CNST_CHAR*,  CNST_CHAR*, SELE2_STYPE , CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      ITYPE*,      STYPE*,      STYPE*,      STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      ITYPE*,      ITYPE*);
        extern VOID		DGEES		( CNST_CHAR*,  CNST_CHAR*, SELE2_DTYPE , CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      ITYPE*,      DTYPE*,      DTYPE*,      DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      ITYPE*,      ITYPE*);
        extern VOID		CGEES		( CNST_CHAR*,  CNST_CHAR*, SELE1_CTYPE , CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      ITYPE*,      CTYPE*,      CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      STYPE*,      ITYPE*,      ITYPE*);
        extern VOID		ZGEES		( CNST_CHAR*,  CNST_CHAR*, SELE1_ZTYPE , CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      ITYPE*,      ZTYPE*,      ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      DTYPE*,      ITYPE*,      ITYPE*);
        // getrf
        extern VOID		SGETRF		(CNST_ITYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      ITYPE*,      ITYPE*);
        extern VOID		DGETRF		(CNST_ITYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      ITYPE*,      ITYPE*);
        extern VOID		CGETRF		(CNST_ITYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      ITYPE*,      ITYPE*);
        extern VOID		ZGETRF		(CNST_ITYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      ITYPE*,      ITYPE*);
        // getri
        extern VOID		SGETRI		(CNST_ITYPE*,      STYPE*, CNST_ITYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID		DGETRI		(CNST_ITYPE*,      DTYPE*, CNST_ITYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID		CGETRI		(CNST_ITYPE*,      CTYPE*, CNST_ITYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID		ZGETRI		(CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      ITYPE*);
        // [sy|he]evd
        extern VOID		SSYEVD		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      STYPE*,      STYPE*, CNST_ITYPE*,      ITYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID		DSYEVD		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      DTYPE*,      DTYPE*, CNST_ITYPE*,      ITYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID		CHEEVD		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      STYPE*,      CTYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      ITYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID		ZHEEVD		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      DTYPE*,      ZTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      ITYPE*, CNST_ITYPE*,      ITYPE*);
        // geev
        extern VOID		SGEEV		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      STYPE*,      STYPE*,      STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID		DGEEV		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      DTYPE*,      DTYPE*,      DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID		CGEEV		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      CTYPE*,      CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      STYPE*,      ITYPE*);
        extern VOID		ZGEEV		( CNST_CHAR*,  CNST_CHAR*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      ZTYPE*,      ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      DTYPE*,      ITYPE*);
        // gesv
        extern VOID     SGESV       (CNST_ITYPE*, CNST_ITYPE*,      STYPE*, CNST_ITYPE*,      ITYPE*,      STYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID     DGESV       (CNST_ITYPE*, CNST_ITYPE*,      DTYPE*, CNST_ITYPE*,      ITYPE*,      DTYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID     CGESV       (CNST_ITYPE*, CNST_ITYPE*,      CTYPE*, CNST_ITYPE*,      ITYPE*,      CTYPE*, CNST_ITYPE*,      ITYPE*);
        extern VOID     ZGESV       (CNST_ITYPE*, CNST_ITYPE*,      ZTYPE*, CNST_ITYPE*,      ITYPE*,      ZTYPE*, CNST_ITYPE*,      ITYPE*);

        // -- LAPACK END-- //
	};
};
#endif // _MATRIXOP_CONFIG_HPP
