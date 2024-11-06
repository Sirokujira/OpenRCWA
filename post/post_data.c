/*
post_data.c

read post data
*/

#include "orcwa.h"
#include "orcwa_prototype.h"

#define MAXTOKEN 1000

int post_data(FILE *fp)
{
	int    ntoken, nline;
	int    version = 0;
	char   prog[BUFSIZ];
	char   strline[BUFSIZ], strkey[BUFSIZ], strsave[BUFSIZ];
	char   *token[MAXTOKEN];
	const char sep[] = " \t";			// separator
	const char errfmt3[] = "*** invalid %s data #%d\n";

	// initialize

	MatchingLoss = 0;

	Piter =
	Pfeed =
	Ppoint = 0;

	IFreq[0] =
	IFreq[1] =
	IFreq[2] =
	IFreq[3] =
	IFreq[4] =
	IFreq[5] = 0;
	Freqdiv = 10;

	IFar0d =
	NFar1d =
	NFar2d =
	NNear1d =
	NNear2d = 0;

	Far1dScale.db = 1;       // dB
	Far1dScale.user = 0;     // auto scale
	Far1dStyle = 0;
	Far1dNorm = 0;
	Far1dComp[0] = 1;
	Far1dComp[1] = 0;
	Far1dComp[2] = 0;

	Far2dScale.db = 1;       // dB
	Far2dScale.user = 0;     // auto scale
	Far2dComp[0] = 1;
	Far2dComp[1] = 0;
	Far2dComp[2] = 0;
	Far2dComp[3] = 0;
	Far2dComp[4] = 0;
	Far2dComp[5] = 0;
	Far2dComp[6] = 0;
	Far2dObj = 0.5;

	Near1dScale.db = 0;      // V/m
	Near1dScale.user = 0;    // auto scale
	Near1dNoinc = 0;

	Near2dDim[0] = Near2dDim[1] = 1;
	Near2dFrame = 0;
	Near2dScale.db = 0;      // V/m
	Near2dScale.user = 0;    // auto scale
	Near2dObj = 1;
	Near2dNoinc = 0;
	Near2dIzoom = 0;

	// window size 2D

	Width2d    = 1.5;
	Height2d   = 1.0;
	Fontsize2d = 0.03;

	// HTML : pixel (for Linux)
	if (HTML) {
		Width2d    = 750;
		Height2d   = 500;
		Fontsize2d = 13;
	}

	// read

	nline = 0;
	while (fgets(strline, sizeof(strline), fp) != NULL) {
		// skip a empty line
		if (strlen(strline) <= 1) continue;

		// skip a comment line
		if (strline[0] == '#') continue;

		// delete "\n"
		//printf("%zd\n", strlen(strline));
		if (strstr(strline, "\r\n") != NULL) {
			strline[strlen(strline) - 2] = '\0';
		}
		else if ((strstr(strline, "\r") != NULL) || (strstr(strline, "\n") != NULL)) {
			strline[strlen(strline) - 1] = '\0';
		}
		//printf("%zd\n", strlen(strline));

		// "end" -> break
		if (!strncmp(strline, "end", 3)) break;

		// save "strline"
		strcpy(strsave, strline);

		// token ("strline" is destroyed)
		ntoken = tokenize(strline, sep, token, MAXTOKEN);
		//for (int i = 0; i < ntoken; i++) printf("%d %s\n", i, token[i]);

		// check number of data and "=" (exclude header)
		if ((nline > 0) && ((ntoken < 3) || strcmp(token[1], "="))) continue;

		// keyword
		strcpy(strkey, token[0]);

		// input
		if      (nline == 0) {
			strcpy(prog, strkey);
			if (strcmp(prog, "OpenRCWA") && strcmp(prog, "OpenTHFD")) {
				printf("%s\n", "*** not OpenRCWA/OpenTHFD data");
				return 1;
			}
			if (ntoken < 3) {
				printf("%s\n", "*** no version data");
				return 1;
			}
			version = (10 * atoi(token[1])) + atoi(token[2]);
			version = version;  // dummy
			nline++;
		}
		else if (!strcmp(strkey, "matchingloss")) {
			MatchingLoss = atoi(token[2]);
		}
		else if (!strcmp(strkey, "plotiter")) {
			Piter = atoi(token[2]);
		}
		else if (!strcmp(strkey, "plotfeed")) {
			Pfeed = atoi(token[2]);
		}
		else if (!strcmp(strkey, "plotpoint")) {
			Ppoint = atoi(token[2]);
		}
		else if (!strcmp(strkey, "plotsmith")) {
			IFreq[0] = atoi(token[2]);
		}
		else if (!strcmp(strkey, "plotzin"     ) ||
		         !strcmp(strkey, "plotyin"     ) ||
		         !strcmp(strkey, "plotref"     ) ||
		         !strcmp(strkey, "plotspara"   ) ||
		         !strcmp(strkey, "plotcoupling")) {
			const int id = !strcmp(strkey, "plotzin"     ) ? 1
			             : !strcmp(strkey, "plotyin"     ) ? 2
			             : !strcmp(strkey, "plotref"     ) ? 3
			             : !strcmp(strkey, "plotspara"   ) ? 4
			             : !strcmp(strkey, "plotcoupling") ? 5 : 0;
			if ((ntoken > 2) && !strcmp(token[2], "1")) {
				IFreq[id] = 1;
				FreqScale[id].user = 0;
			}
			else if ((ntoken > 5) && !strcmp(token[2], "2")) {
				IFreq[id] = 1;
				FreqScale[id].user = 1;
				FreqScale[id].min = atof(token[3]);
				FreqScale[id].max = atof(token[4]);
				FreqScale[id].div = atoi(token[5]);
			}
		}
		else if (!strcmp(strkey, "freqdiv")) {
			Freqdiv = atoi(token[2]);
		}
		else if (!strcmp(strkey, "plotfar0d")) {
			if ((ntoken > 4) && !strcmp(token[4], "1")) {
				IFar0d = 1;
				Far0d[0] = atof(token[2]);
				Far0d[1] = atof(token[3]);
				Far0dScale.user = 0;
			}
			else if ((ntoken > 7) && !strcmp(token[4], "2")) {
				IFar0d = 1;
				Far0d[0] = atof(token[2]);
				Far0d[1] = atof(token[3]);
				Far0dScale.user = 1;
				Far0dScale.min = atof(token[5]);
				Far0dScale.max = atof(token[6]);
				Far0dScale.div = atoi(token[7]);
			}
		}
		else if (!strcmp(strkey, "plotfar1d")) {
			char dir = (char)toupper((int)token[2][0]);
			if ((dir != 'X') && (dir != 'Y') && (dir != 'Z') &&
			    (dir != 'V') && (dir != 'H')) {
				printf(errfmt3, strkey, NFar1d + 1);
				return 1;
			}
			if ((((dir == 'X') || (dir == 'Y') || (dir == 'Z')) && (ntoken < 4)) ||
			    (((dir == 'V') || (dir == 'H')) && (ntoken < 5))) {
				printf(errfmt3, strkey, NFar1d + 1);
				return 1;
			}
			Far1d = (far1d_t *)realloc(Far1d, (NFar1d + 1) * sizeof(far1d_t));
			Far1d[NFar1d].dir = dir;
			Far1d[NFar1d].div = atoi(token[3]);
			if ((dir == 'V') || (dir == 'H')) {
				Far1d[NFar1d].angle = atof(token[4]);
			}
			NFar1d++;
		}
		else if (!strcmp(strkey, "plotfar2d")) {
			if (ntoken > 3) {
				Far2d.divtheta = atoi(token[2]);
				Far2d.divphi   = atoi(token[3]);
				NFar2d = 1;
			}
		}
		else if (!strcmp(strkey, "plotnear1d")) {
			if ((ntoken < 6) || (strlen(token[2]) > 2) || (strlen(token[3]) > 1)) {
				printf(errfmt3, strkey, NNear1d + 1);
				return 1;
			}
			Near1d = (near1d_t *)realloc(Near1d, (NNear1d + 1) * sizeof(near1d_t));
			strcpy(Near1d[NNear1d].cmp, token[2]);
			Near1d[NNear1d].dir = (char)toupper((int)token[3][0]);
			Near1d[NNear1d].pos1 = atof(token[4]);
			Near1d[NNear1d].pos2 = atof(token[5]);
			NNear1d++;
		}
		else if (!strcmp(strkey, "plotnear2d")) {
			if ((ntoken < 5) || (strlen(token[2]) > 2) || (strlen(token[3]) > 1)) {
				printf(errfmt3, strkey, NNear2d + 1);
				return 1;
			}
			Near2d = (near2d_t *)realloc(Near2d, (NNear2d + 1) * sizeof(near2d_t));
			//pos2d0 = (double *)realloc(pos2d0, (NNear2d + 1) * sizeof(double));
			strcpy(Near2d[NNear2d].cmp, token[2]);
			Near2d[NNear2d].dir = (char)toupper((int)token[3][0]);
			Near2d[NNear2d].pos0 = atof(token[4]);
			NNear2d++;
		}
		else if (!strcmp(strkey, "far1dcomponent")) {
			if (ntoken > 4) {
				for (int n = 0; n < 3; n++) {
					Far1dComp[n] = atoi(token[2 + n]);
				}
			}
		}
		else if (!strcmp(strkey, "far1dstyle")) {
			Far1dStyle = atoi(token[2]);
		}
		else if (!strcmp(strkey, "far1ddb")) {
			Far1dScale.db = atoi(token[2]);
		}
		else if (!strcmp(strkey, "far1dnorm")) {
			Far1dNorm = atoi(token[2]);
		}
		else if (!strcmp(strkey, "far1dscale")) {
			if (ntoken > 4) {
				Far1dScale.user = 1;
				Far1dScale.min = atof(token[2]);
				Far1dScale.max = atof(token[3]);
				Far1dScale.div = atoi(token[4]);
			}
		}
		else if (!strcmp(strkey, "far2dcomponent")) {
			if (ntoken > 8) {
				for (int n = 0; n < 7; n++) {
					Far2dComp[n] = atoi(token[2 + n]);
				}
			}
		}
		else if (!strcmp(strkey, "far2ddb")) {
			Far2dScale.db = atoi(token[2]);
		}
		else if (!strcmp(strkey, "far2dscale")) {
			if (ntoken > 3) {
				Far2dScale.user = 1;
				Far2dScale.min  = atof(token[2]);
				Far2dScale.max  = atof(token[3]);
			}
		}
		else if (!strcmp(strkey, "far2dobj")) {
			Far2dObj = atof(token[2]);
		}
		else if (!strcmp(strkey, "near1ddb")) {
			Near1dScale.db = atoi(token[2]);
		}
		else if (!strcmp(strkey, "near1dscale")) {
			if (ntoken > 4) {
				Near1dScale.user = 1;
				Near1dScale.min = atof(token[2]);
				Near1dScale.max = atof(token[3]);
				Near1dScale.div = atoi(token[4]);
			}
		}
		else if (!strcmp(strkey, "near1dnoinc")) {
			Near1dNoinc = atoi(token[2]);
		}
		else if (!strcmp(strkey, "near2ddim")) {
			if (ntoken > 3) {
				Near2dDim[0] = atoi(token[2]);
				Near2dDim[1] = atoi(token[3]);
			}
		}
		else if (!strcmp(strkey, "near2dframe")) {
			Near2dFrame = atoi(token[2]);
		}
		else if (!strcmp(strkey, "near2ddb")) {
			Near2dScale.db = atoi(token[2]);
		}
		else if (!strcmp(strkey, "near2dscale")) {
			if (ntoken > 3) {
				Near2dScale.user = 1;
				Near2dScale.min = atof(token[2]);
				Near2dScale.max = atof(token[3]);
				// token[4] : not used
			}
		}
		else if (!strcmp(strkey, "near2dcontour")) {
			Near2dContour = atoi(token[2]);
		}
		else if (!strcmp(strkey, "near2dobj")) {
			Near2dObj = atoi(token[2]);
		}
		else if (!strcmp(strkey, "near2dnoinc")) {
			Near2dNoinc = atoi(token[2]);
		}
		else if (!strcmp(strkey, "near2dzoom")) {
			if (ntoken > 5) {
				Near2dIzoom = 1;
				Near2dHzoom[0] = MIN(atof(token[2]), atof(token[3]));
				Near2dHzoom[1] = MAX(atof(token[2]), atof(token[3]));
				Near2dVzoom[0] = MIN(atof(token[4]), atof(token[5]));
				Near2dVzoom[1] = MAX(atof(token[4]), atof(token[5]));
			}
		}
	}

	// error check

	for (int n = 0; n < NFar1d; n++) {
		if (Far1d[n].div < 2) {
			printf("%s #%d\n", "*** invalid far1d division", n + 1);
			return 1;
		}
	}
	if (NFar2d) {
		if ((Far2d.divtheta <= 0) || (Far2d.divphi <= 0)) {
			printf("%s\n", "*** invalid far2d division");
			return 1;
		}
	}
	for (int n = 0; n < NNear1d; n++) {
		if ((Near1d[n].dir != 'X') && (Near1d[n].dir != 'Y') && (Near1d[n].dir != 'Z')) {
			printf("%s #%d\n", "*** invalid near1d direction", n + 1);
			return 1;
		}
	}
	for (int n = 0; n < NNear2d; n++) {
		if ((Near2d[n].dir != 'X') && (Near2d[n].dir != 'Y') && (Near2d[n].dir != 'Z')) {
			printf("%s #%d\n", "*** invalid near2d direction", n + 1);
			return 1;
		}
	}

	return 0;
}
