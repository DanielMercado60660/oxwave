/* Copyright (c) 2008, Colorado School of Mines, */
/* All rights reserved. */
/* */
/* Redistribution and use in source and binary forms, with or */
/* without modification, are permitted provided that the following */
/* conditions are met: */
/* */
/*     *  Redistributions of source code must retain the above copyright */
/*        notice, this list of conditions and the following disclaimer. */
/*     *  Redistributions in binary form must reproduce the above */
/*        copyright notice, this list of conditions and the following */
/*        disclaimer in the documentation and/or other materials provided */
/*        with the distribution. */
/*     *  Neither the name of the Colorado School of Mines nor the names of */
/*        its contributors may be used to endorse or promote products */
/*        derived from this software without specific prior written permission. */
/* */
/* Warranty Disclaimer: */
/* THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE */
/* COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, */
/* STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING */
/* IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */
/* */
/* Rust port for OxWave, a Rust port of CWP/SU. */

/*
 * Test shim for oxwave-cref: copies of functions that are `static` inside
 * SU programs (so not in any .a), made callable for differential tests.
 * ibm_to_float: verbatim from su/main/data_conversion/segyread.c
 * float_to_ibm: verbatim from su/main/data_conversion/segywrite.c (SUXDR branch)
 * Compiled with the reference SU build's flags (-O0 -ffp-contract=off).
 */
#include <string.h>
#include <unistd.h>

/* libpar's pagedoc() references the program's sdoc[] */
char *sdoc[] = { "oxwave-cref test harness", 0 };

typedef int fourbyte;

void oxcref_ibm_to_float(int from[], int to[], int n, int endian, int verbose)
{
    /* Initialized data */

    static int first = 1;

    /* System generated locals */
    int i__1;
    int j,k;

    /* Local variables */
    int   *in;
    float *out;
    int eibm, i__, mhibm;
    static int m1[512];
    static float r1[512];

    unsigned int jj;
    union {
      float rrf;
      int iif;
      unsigned int uuf;
    } cvtmp;

    float r_infinity__;
    int et3e;

    (void) verbose;
    if(endian == 0) {
      swab(from,to,n*sizeof(int));
      for(i__ = 0; i__<n; ++i__) {
        j = to[i__];
        k = j<<16;
        to[i__] = k+((j>>16)&65535);
      }
      in = to;
    } else {
      in = from;
    }
    /* Parameter adjustments */
    out = (float *) to;
    --out;
    --in;
    /* Function Body */

    if (first) {
	first = ! first;
	cvtmp.iif = 2139095039;
	r_infinity__ = cvtmp.rrf;
	for (i__ = 0; i__ <= 511; ++i__) {
	    i__1 = i__ & 255;
	    eibm = i__1 >> 1;
	    mhibm = i__ & 1;
	    et3e = (eibm << 2) - 130;
	    if (et3e > 0 && et3e <= 255) {
		i__1 = et3e ^ (i__ & 255);
		m1[i__] = i__1 << 23;
		if (mhibm == 1) {
		    r1[i__] = 0.f;
		} else {
		    i__1 = et3e | (i__ & 256);
		    cvtmp.iif = i__1 << 23;
		    r1[i__] = -(cvtmp.rrf);
		}
	    } else if (et3e <= 0) {
		m1[i__] = i__ << 23;
		r1[i__] = 0.f;
	    } else {
		m1[i__] = i__ << 23;
		if (i__ < 256) {
		    r1[i__] = r_infinity__;
		} else {
		    r1[i__] = -r_infinity__;
		}
	    }
	}
    }

    for (i__ = 1; i__ <= n; ++i__) {
	cvtmp.iif = in[i__];
	jj = cvtmp.uuf>>23;
	cvtmp.iif = m1[jj] ^ cvtmp.iif;
	out[i__] = cvtmp.rrf + r1[jj];
    }
}

void oxcref_float_to_ibm(fourbyte *from, fourbyte *to, int n, int endian)
{
    register fourbyte fconv, fmant, t;
    register int i;

    for (i=0;i<n;++i) {
	fconv = from[i];
	if (fconv) {
	    fmant = (0x007fffff & fconv) | 0x00800000;
	    t = (fourbyte) ((0x7f800000 & fconv) >> 23) - 126;
	    while (t & 0x3) { ++t; fmant >>= 1; }
	    fconv = (0x80000000 & fconv) | (((t>>2) + 64) << 24) | fmant;
	}
	if(endian==0)
		fconv = (fconv<<24) | ((fconv>>24)&0xff) |
			((fconv&0xff00)<<8) | ((fconv&0xff0000)>>8);

	to[i] = fconv;
    }
    return;
}
