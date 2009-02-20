/* FUNPACK
 * R. Seaman, NOAO
 * uses fits_img_compress by W. Pence, HEASARC
 */

#include "fitsio.h"
#include "fpack.h"

int main (int argc, char *argv[])
{
	fpstate	fpvar;

	if (argc <= 1) { fu_usage (); fu_hint (); exit (-1); }

	fp_init (&fpvar);
	fu_get_param (argc, argv, &fpvar);

	if (fpvar.listonly) {
	    fp_list (argc, argv, fpvar);

	} else {
	    fp_preflight (argc, argv, &fpvar);
	    fp_loop (argc, argv, FUNPACK, fpvar);
	}

	exit (0);
}

fu_get_param (int argc, char *argv[], fpstate *fpptr)
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

		if (argv[iarg][1] == 'F') {
		    fpptr->clobber++;

		} else if (argv[iarg][1] == 'A') {
		    if (++iarg >= argc) {
			fu_usage (); fu_hint (); exit (-1);
		    } else
			strncpy (fpptr->suffix, argv[iarg], SZ_STR);

		} else if (argv[iarg][1] == 'P') {
		    if (++iarg >= argc) {
			fu_usage (); fu_hint (); exit (-1);
		    } else
			strncpy (fpptr->prefix, argv[iarg], SZ_STR);

		} else if (argv[iarg][1] == 'S') {
		    fpptr->to_stdout++;

		} else if (argv[iarg][1] == 'L') {
		    fpptr->listonly++;

		} else if (argv[iarg][1] == 'C') {
		    fpptr->do_checksums = 0;

		} else if (argv[iarg][1] == 'H') {
		    fu_help (); exit (0);

		} else if (argv[iarg][1] == 'V') {
		    fp_version (); exit (0);

		} else {
		    fp_msg ("unknown command line flag `");
		    fp_msg (argv[iarg]); fp_msg ("'\n");
		    fu_usage (); fu_hint (); exit (-1);
		}

	    } else
		break;
	}

	if (! fpptr->listonly && ! fpptr->to_stdout && ! fpptr->clobber &&
	    ! fpptr->prefix[0] && ! fpptr->suffix[0]) {

	    fp_msg ("to overwrite input files, must specify `-F'\n");
            fp_msg ("  otherwise specify output suffix with `-A'\n");
            fp_msg ("  or output prefix with `-P'\n\n");

	    fu_usage (); fu_hint (); exit (-1);
	}

	if (iarg >= argc) {
	    fp_msg ("no FITS files to uncompress\n");
	    fu_usage (); exit (-1);
	} else
	    fpptr->firstfile = iarg;
}

fu_usage ()
{
	fp_msg ("usage: funpack [-F] [-P <pre>|-A <suffix>] [-S] [-L] [-C] [-H] [-V] <FITS>\n");
}

fu_hint ()
{
	fp_msg ("      `funpack -H' for help\n");
}

fu_help ()
{
fp_msg ("funpack, decompress fpacked files.  Version ");
fp_version ();
fu_usage ();
fp_msg ("\n");

fp_msg ("Flags must be separate and appear before filenames:\n");
fp_msg ("   -v          verbose list of each file\n");
fp_msg ("   -F          clobber output [required to overwrite in-place]\n");
fp_msg ("   -P <pre>    prepend <pre> to create new output filenames\n");
fp_msg ("   -A <suffix> append <suffix> to create new output filenames\n");
fp_msg ("   -S          output to STDOUT\n");
fp_msg ("   -L          list contents, files unchanged\n");

fp_msg ("   -C          don't update FITS checksum keywords\n");

fp_msg ("   -H          print this message\n");
fp_msg ("   -V          print version number\n");

fp_msg (" <FITS>        FITS files to unpack\n");
}
