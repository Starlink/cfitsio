/* used by FPACK and FUNPACK
 * R. Seaman, NOAO
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#define	FPACK_VERSION	"0.9.6 (April 2008)"
#define	FP_INIT_MAGIC	42
#define	FPACK		0
#define	FUNPACK		1

#define	DEF_QLEVEL	16.
#define	DEF_HCOMP_SCALE	 0.
#define	DEF_HCOMP_SMOOTH 0
#define	DEF_RESCALE_NOISE 0

#define	SZ_STR		161
#define	SZ_CARD		81

typedef struct
{
	int	comptype;
	float	quantize_level;
	float	scale;
	float	rescale_noise;
	int	smooth;
	long	ntile[MAX_COMPRESS_DIM];
	int	islossless;

	int	to_stdout;
	int	listonly;
	int	clobber;
	int	do_checksums;
	int	test_all;
	int	verbose;

	char	prefix[SZ_STR];
	char	suffix[SZ_STR];
	char	outfile[SZ_STR];
	int	firstfile;

	int	initialized;
	int	preflight_checked;
} fpstate;

typedef struct
{
	int	n_nulls;
	double	minval;
	double 	maxval;
	double	mean;
	double	sigma;
	double	noise1;
	double	noise3;
} imgstats;
