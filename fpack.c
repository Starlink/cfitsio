/* FPACK
 * R. Seaman, NOAO, with a few enhancements by W. Pence, HEASARC
 *
 * Calls fits_img_compress in the CFITSIO library by W. Pence, HEASARC
 */

#include "fitsio.h"
#include "fpack.h"

int main(int argc, char *argv[])
{
	fpstate	fpvar;

	if (argc <= 1) { fp_usage (); fp_hint (); exit (-1); }

	fp_init (&fpvar);
	fp_get_param (argc, argv, &fpvar);

	if (fpvar.listonly) {
	    fp_list (argc, argv, fpvar);

	} else {
	    fp_preflight (argc, argv, &fpvar);
	    fp_loop (argc, argv, FPACK, fpvar);
	}

	exit (0);
}

fp_get_param (int argc, char *argv[], fpstate *fpptr)
{
	int	gottype=0, gottile=0, wholetile=0, iarg, len, ndim, ii;
	char	tmp[SZ_STR], tile[SZ_STR];

        if (fpptr->initialized != FP_INIT_MAGIC) {
            fp_msg ("internal initialization error\n"); exit (-1);
        }

	tile[0] = (char) NULL;

	/* flags must come first and be separately specified
	 */
	for (iarg = 1; iarg < argc; iarg++) {
	    if (argv[iarg][0] == '-' && strlen (argv[iarg]) == 2) {

		/* Rice is the default, so -r is superfluous 
		 */
		if (argv[iarg][1] == 'r') {
		    fpptr->comptype = RICE_1;
		    if (gottype) {
			fp_msg ("multiple compression flags\n");
			fp_usage (); exit (-1);
		    } else
			gottype++;

		} else if (argv[iarg][1] == 'p') {
		    fpptr->comptype = PLIO_1;
		    if (gottype) {
			fp_msg ("multiple compression flags\n");
			fp_usage (); exit (-1);
		    } else
			gottype++;

		} else if (argv[iarg][1] == 'g') {
		    fpptr->comptype = GZIP_1;
		    if (gottype) {
			fp_msg ("multiple compression flags\n");
			fp_usage (); exit (-1);
		    } else
			gottype++;

		} else if (argv[iarg][1] == 'h') {
		    fpptr->comptype = HCOMPRESS_1;
		    if (gottype) {
			fp_msg ("multiple compression flags\n");
			fp_usage (); exit (-1);
		    } else
			gottype++;

		} else if (argv[iarg][1] == 'd') {
		    fpptr->comptype = NOCOMPRESS;
		    if (gottype) {
			fp_msg ("multiple compression flags\n");
			fp_usage (); exit (-1);
		    } else
			gottype++;

		} else if (argv[iarg][1] == 'q') {
		    if (++iarg >= argc) {
			fp_usage (); exit (-1);
		    } else {
			fpptr->quantize_level = (float) atof (argv[iarg]);
		    }
		} else if (argv[iarg][1] == 'n') {
		    if (++iarg >= argc) {
			fp_usage (); exit (-1);
		    } else {
			fpptr->rescale_noise = (float) atof (argv[iarg]);
		    }
		} else if (argv[iarg][1] == 's') {
		    if (++iarg >= argc) {
			fp_usage (); exit (-1);
		    } else {
			fpptr->scale = (float) atof (argv[iarg]);
		    }
		} else if (argv[iarg][1] == 't') {
		    if (gottile) {
			fp_msg ("multiple tile specifications\n");
			fp_usage (); exit (-1);
		    } else
			gottile++;

		    if (++iarg >= argc) {
			fp_usage (); exit (-1);
		    } else
			strncpy (tile, argv[iarg], SZ_STR); /* checked below */

		} else if (argv[iarg][1] == 'v') {
		    fpptr->verbose = 1;

		} else if (argv[iarg][1] == 'w') {
		    wholetile++;
		    if (gottile) {
			fp_msg ("multiple tile specifications\n");
			fp_usage (); exit (-1);
		    } else
			gottile++;

		} else if (argv[iarg][1] == 'F') {
		    fpptr->clobber++;

		} else if (argv[iarg][1] == 'A') {
		    if (++iarg >= argc) {
			fp_usage (); fp_hint (); exit (-1);
		    } else
			strncpy (fpptr->suffix, argv[iarg], SZ_STR);

		} else if (argv[iarg][1] == 'P') {
		    if (++iarg >= argc) {
			fp_usage (); fp_hint (); exit (-1);
		    } else
			strncpy (fpptr->prefix, argv[iarg], SZ_STR);

		} else if (argv[iarg][1] == 'O') {
		    if (++iarg >= argc) {
			fp_usage (); fp_hint (); exit (-1);
		    } else
			strncpy (fpptr->outfile, argv[iarg], SZ_STR);

		} else if (argv[iarg][1] == 'S') {
		    fpptr->to_stdout++;

		} else if (argv[iarg][1] == 'L') {
		    fpptr->listonly++;

		} else if (argv[iarg][1] == 'C') {
		    fpptr->do_checksums = 0;

		} else if (argv[iarg][1] == 'T') {
		    fpptr->test_all = 1;

		} else if (argv[iarg][1] == 'H') {
		    fp_help (); exit (0);

		} else if (argv[iarg][1] == 'V') {
		    fp_version (); exit (0);

		} else {
		    fp_msg ("unknown command line flag `");
		    fp_msg (argv[iarg]); fp_msg ("'\n");
		    fp_usage (); fp_hint (); exit (-1);
		}

	    } else
		break;
	}

        if (! fpptr->listonly && ! fpptr->to_stdout && ! fpptr->clobber &&
	    ! fpptr->prefix[0] && ! fpptr->suffix[0] && ! fpptr->test_all) {

            fp_msg ("to overwrite input files, must specify `-F'\n");
            fp_msg ("  otherwise specify output suffix with `-A'\n");
            fp_msg ("  or output prefix with `-P'\n\n");

            fp_usage (); fp_hint (); exit (-1);
        }


	if (fpptr->scale != 0. && 
	         fpptr->comptype != HCOMPRESS_1 && fpptr->test_all != 1) {

	    fp_msg ("`-s' requires `-h or -T'\n"); exit (-1);
	}


        if (fpptr->clobber && fpptr->scale != 0. &&
		! fpptr->prefix[0] && ! fpptr->suffix[0] ) {

            fp_msg ("Can only overwrite input files when performing\n");
            fp_msg ("  lossless compression  (s must equal 0.)\n\n");

            fp_usage (); fp_hint (); exit (-1);
        }

	if (wholetile) {
	    for (ndim=0; ndim < MAX_COMPRESS_DIM; ndim++)
		fpptr->ntile[ndim] = (long) 0;

	} else if (gottile) {
	    len = strlen (tile);
	    for (ii=0, ndim=0; ii < len; ) {
		if (! (isdigit (tile[ii]) || tile[ii] == ',')) {
		    fp_msg ("`-t' requires comma separated tile dims, ");
		    fp_msg ("e.g., `-t 100,100'\n"); exit (-1);
		}

		if (tile[ii] == ',') { ii++; continue; }

		fpptr->ntile[ndim] = atol (&tile[ii]);
		for ( ; isdigit(tile[ii]); ii++);

		if (++ndim > MAX_COMPRESS_DIM) {
		    fp_msg ("too many dimensions for `-t', max=");
		    sprintf (tmp, "%d\n", MAX_COMPRESS_DIM); fp_msg (tmp);
		    exit (-1);
		}
	    }
	}

	if (iarg >= argc) {
	    fp_msg ("no FITS files to compress\n");
	    fp_usage (); exit (-1);
	} else
	    fpptr->firstfile = iarg;
}

fp_usage ()
{
fp_msg ("usage: fpack ");
fp_msg (
"[-r|-h|-g|-p] [-w|-t <axes>] [-q <level>] [-s <scale>] [-n <noise>] -v <FITS>\n");
fp_msg ("more:   [-T] [-F] [-P <pre>|-A <suffix>] [-S] [-L] [-C] [-H] [-V]\n");
}

fp_hint () { fp_msg ("      `fpack -H' for help\n"); }

fp_help ()
{
fp_msg ("fpack, a FITS tile-compression engine.  Version ");
fp_version ();
fp_usage ();
fp_msg ("\n");

fp_msg ("Flags must be separate and appear before filenames:\n");
fp_msg ("   -r          Rice compression [default], or\n");
fp_msg ("   -h          Hcompress compression, or\n");
fp_msg ("   -g          GZIP (per-tile) compression, or\n");
fp_msg ("   -p          PLIO compression (only for positive 8 or 16-bit integer images)\n");
fp_msg ("   -d          tile the image without compression (debugging mode)\n");

fp_msg ("   -w          compress the whole image,as a single large tile\n");
fp_msg ("   -t <axes>   comma separated list of tile sizes [default=row]\n");
fp_msg ("   -q <level>  quantization level for real pixels [default=16]\n");
fp_msg ("               (+values relative to RMS noise; -value is absolute)\n");

fp_msg ("   -s <scale>  scale factor for lossy Hcompress [default = 0 = lossless]\n");
fp_msg ("               (+values relative to RMS noise; -value is absolute)\n");
fp_msg ("   -n <noise>  rescale scaled-integer images to reduce noise\n");

fp_msg ("   -v          verbose mode; list each file as it is processed\n");
fp_msg ("   -T          print test comparison report of compression algorithms\n");
fp_msg ("   -O <file>   write output text file with test results\n");

fp_msg ("\nkeywords shared with funpack:\n");

fp_msg ("   -F          clobber output [required to overwrite in-place]\n");
fp_msg ("   -P <pre>    prepend <pre> to create new output filenames\n");
fp_msg ("   -A <suffix> append <suffix> to create new output filenames\n");
fp_msg ("   -S          output compressed FITS files to STDOUT\n");
fp_msg ("   -L          list contents, files unchanged\n");

fp_msg ("   -C          don't update FITS checksum keywords\n");

fp_msg ("   -H          print this message\n");
fp_msg ("   -V          print version number\n");

fp_msg (" <FITS>        FITS files to pack\n");
}
