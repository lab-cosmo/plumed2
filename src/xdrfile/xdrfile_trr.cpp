/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2009-2014, Erik Lindahl & David van der Spoel
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id$
 *
 /* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*-
 *
 * $Id$
 *
 * Copyright (c) 2009-2014, Erik Lindahl & David van der Spoel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "xdrfile.h"
#include "xdrfile_trr.h"

#define BUFSIZE		128
#define GROMACS_MAGIC   1993

namespace PLMD {
namespace xdrfile {

typedef struct		/* This struct describes the order and the	*/
/* sizes of the structs in a trjfile, sizes are given in bytes.	*/
{
    mybool  bDouble;        /* Double precision?                            */
    int	ir_size;	/* Backward compatibility		        */
    int	e_size;		/* Backward compatibility		        */
    int	box_size;	/* Non zero if a box is present			*/
    int   vir_size;       /* Backward compatibility		        */
    int   pres_size;      /* Backward compatibility		        */
    int	top_size;	/* Backward compatibility		        */
    int	sym_size;	/* Backward compatibility		        */
    int	x_size;		/* Non zero if coordinates are present		*/
    int	v_size;		/* Non zero if velocities are present		*/
    int	f_size;		/* Non zero if forces are present		*/

    int	natoms;		/* The total number of atoms			*/
    int	step;		/* Current step number				*/
    int	nre;		/* Backward compatibility		        */
    float	tf;		/* Current time					*/
    float	lambdaf;		/* Current value of lambda			*/
    double	td;		/* Current time					*/
    double	lambdad;		/* Current value of lambda			*/
} t_trnheader;

static int nFloatSize(t_trnheader *sh,int *nflsz)
{
    int nflsize=0;
  
    if (sh->box_size)
        nflsize = sh->box_size/(DIM*DIM);
    else if (sh->x_size)
        nflsize = sh->x_size/(sh->natoms*DIM);
    else if (sh->v_size)
        nflsize = sh->v_size/(sh->natoms*DIM);
    else if (sh->f_size)
        nflsize = sh->f_size/(sh->natoms*DIM);
    else 
        return exdrHEADER;
  
    if (((nflsize != sizeof(float)) && (nflsize != sizeof(double))))
        return exdrHEADER;
      
    *nflsz = nflsize;
  
    return exdrOK;
}

static int do_trnheader(XDRFILE *xd,mybool bRead,t_trnheader *sh)
{
	int magic=GROMACS_MAGIC;
	int nflsz,slen,result;
	const char *version = "GMX_trn_file";
	char buf[BUFSIZE];
  
	if (xdrfile_read_int(&magic,1,xd) != 1)
		return exdrINT;
  
	if (bRead) 
    {
        if (xdrfile_read_int(&slen,1,xd) != 1)
            return exdrINT;
        if (slen != strlen(version)+1)
            return exdrSTRING;
        if (xdrfile_read_string(buf,BUFSIZE,xd) <= 0)
            return exdrSTRING;
    }
	else 
    {
        slen = strlen(version)+1;
        if (xdrfile_read_int(&slen,1,xd) != 1)
            return exdrINT;
        if (xdrfile_write_string(const_cast<char*>(version),xd) != (strlen(version)+1) )
            return exdrSTRING;
    }
	if (xdrfile_read_int(&sh->ir_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->e_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->box_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->vir_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->pres_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->top_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->sym_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->x_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->v_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->f_size,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->natoms,1,xd) != 1)
		return exdrINT;
	
	if ((result = nFloatSize(sh,&nflsz)) != exdrOK)
		return result;
	sh->bDouble = (nflsz == sizeof(double));
	
	if (xdrfile_read_int(&sh->step,1,xd) != 1)
		return exdrINT;
	if (xdrfile_read_int(&sh->nre,1,xd) != 1)
		return exdrINT;
	if (sh->bDouble) 
    {
        if (xdrfile_read_double(&sh->td,1,xd) != 1)
            return exdrDOUBLE;
        sh->tf = sh->td;
        if (xdrfile_read_double(&sh->lambdad,1,xd) != 1)
            return exdrDOUBLE;
        sh->lambdaf = sh->lambdad;
    }
	else 
    {
        if (xdrfile_read_float(&sh->tf,1,xd) != 1)
            return exdrFLOAT;
        sh->td = sh->tf;
        if (xdrfile_read_float(&sh->lambdaf,1,xd) != 1)
            return exdrFLOAT;
        sh->lambdad = sh->lambdaf;
    }
  
    return exdrOK;
}

static int do_htrn(XDRFILE *xd,mybool bRead,t_trnheader *sh,
				   matrix box,rvec *x,rvec *v,rvec *f)
{
	double pvd[DIM*DIM];
	double *dx=NULL;
	float  pvf[DIM*DIM];
	float  *fx=NULL;
	int    i,j;
	
	if (sh->bDouble) 
	{
		if (sh->box_size != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<DIM); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != box)
                        {
                            pvd[i*DIM+j] = box[i][j];
                        }
            }
            if (xdrfile_read_double(pvd,DIM*DIM,xd) == DIM*DIM) 
            {
                for(i=0; (i<DIM); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != box)
                        {
                            box[i][j] = pvd[i*DIM+j];
                        }
            }
            else
                return exdrDOUBLE;
        }
			
		if (sh->vir_size != 0) 
        {
            if (xdrfile_read_double(pvd,DIM*DIM,xd) != DIM*DIM) 
                return exdrDOUBLE;
        }
		
		if (sh->pres_size!= 0) 
        {
            if (xdrfile_read_double(pvd,DIM*DIM,xd) != DIM*DIM) 
                return exdrDOUBLE;
        }
		
		if ((sh->x_size != 0) || (sh->v_size != 0) || (sh->f_size != 0)) {
			dx = (double *)calloc(sh->natoms*DIM,sizeof(dx[0]));
			if (NULL == dx)
				return exdrNOMEM;
		}
		if (sh->x_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            dx[i*DIM+j] = x[i][j];
                        }
            }
            if (xdrfile_read_double(dx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                if (bRead) 
                {
                    for(i=0; (i<sh->natoms); i++)
                        for(j=0; (j<DIM); j++)
                            if (NULL != x)
                            {
                                x[i][j] = dx[i*DIM+j];
                            }
                }
            }
            else
                return exdrDOUBLE;
        }
		if (sh->v_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            dx[i*DIM+j] = v[i][j];
                        }
            }
            if (xdrfile_read_double(dx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != v)
                        {
                            v[i][j] = dx[i*DIM+j];
                        }
            }
            else
                return exdrDOUBLE;
        }
		if (sh->f_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            dx[i*DIM+j] = f[i][j];
                        }
            }
            if (xdrfile_read_double(dx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                for(i=0; (i<sh->natoms); i++)
                {
                    for(j=0; (j<DIM); j++)
                    {
                        if (NULL != f)
                        {
                            f[i][j] = dx[i*DIM+j];
                        }
                    }
                }
            }
            else
                return exdrDOUBLE;
        }
		if ((sh->x_size != 0) || (sh->v_size != 0) || (sh->f_size != 0)) {
			free(dx);
		}
	}
	else
		/* Float */
	{
		if (sh->box_size != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<DIM); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != box)
                        {
                            pvf[i*DIM+j] = box[i][j];
                        }
            }
            if (xdrfile_read_float(pvf,DIM*DIM,xd) == DIM*DIM) 
            {
                for(i=0; (i<DIM); i++)
                {
                    for(j=0; (j<DIM); j++)
                    {
                        if (NULL != box)
                        {
                            box[i][j] = pvf[i*DIM+j];
                        }
                    }
                }
            }
            else
                return exdrFLOAT;
        }
			
		if (sh->vir_size != 0) 
        {
            if (xdrfile_read_float(pvf,DIM*DIM,xd) != DIM*DIM) 
                return exdrFLOAT;
        }
		
		if (sh->pres_size!= 0) 
        {
            if (xdrfile_read_float(pvf,DIM*DIM,xd) != DIM*DIM) 
                return exdrFLOAT;
        }
		
		if ((sh->x_size != 0) || (sh->v_size != 0) || (sh->f_size != 0)) {
			fx = (float *)calloc(sh->natoms*DIM,sizeof(fx[0]));
			if (NULL == fx)
				return exdrNOMEM;
		}
		if (sh->x_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            fx[i*DIM+j] = x[i][j];
                        }
            }
            if (xdrfile_read_float(fx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                if (bRead) 
                {
                    for(i=0; (i<sh->natoms); i++)
                        for(j=0; (j<DIM); j++)
                            if (NULL != x)
                                x[i][j] = fx[i*DIM+j];
                }
            }
            else
                return exdrFLOAT;
        }
		if (sh->v_size   != 0) 
        {
            if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            fx[i*DIM+j] = v[i][j];
                        }
            }
            if (xdrfile_read_float(fx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != v)
                            v[i][j] = fx[i*DIM+j];
            }
            else
                return exdrFLOAT;
        }
		if (sh->f_size   != 0) 
        {
           if (!bRead) 
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != x)
                        {
                            fx[i*DIM+j] = f[i][j];
                        }
            }
             if (xdrfile_read_float(fx,sh->natoms*DIM,xd) == sh->natoms*DIM)
            {
                for(i=0; (i<sh->natoms); i++)
                    for(j=0; (j<DIM); j++)
                        if (NULL != f)
                            f[i][j] = fx[i*DIM+j];
            }
            else
                return exdrFLOAT;
        }
		if ((sh->x_size != 0) || (sh->v_size != 0) || (sh->f_size != 0)) {
			free(fx);
		}
	}
	return exdrOK;
}

static int do_trn(XDRFILE *xd,mybool bRead,int *step,float *t,float *lambda,
				  matrix box,int *natoms,rvec *x,rvec *v,rvec *f)
{
    t_trnheader *sh;
    int result;
  
    sh = (t_trnheader *)calloc(1,sizeof(*sh));
  
    if (!bRead) {
        sh->box_size = (NULL != box) ? sizeof(matrix):0;
        sh->x_size   = ((NULL != x) ? (*natoms*sizeof(x[0])):0);
        sh->v_size   = ((NULL != v) ? (*natoms*sizeof(v[0])):0);
        sh->f_size   = ((NULL != f) ? (*natoms*sizeof(f[0])):0);
        sh->natoms = *natoms;
        sh->step   = *step;
        sh->nre    = 0;
        sh->td      = *t;
        sh->lambdad = *lambda;
        sh->tf      = *t;
        sh->lambdaf = *lambda;
    }
    if ((result = do_trnheader(xd,bRead,sh)) != exdrOK)
        return result;
    if (bRead) {
        *natoms = sh->natoms;
        *step   = sh->step;
        *t      = sh->td;
        *lambda = sh->lambdad;
    }
    if ((result = do_htrn(xd,bRead,sh,box,x,v,f)) != exdrOK)
        return result;

    free(sh);
  
    return exdrOK;
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/
 
int read_trr_natoms(char *fn,int *natoms)
{
	XDRFILE *xd;
	t_trnheader sh;
	int  result;
	
	xd = xdrfile_open(fn,"r");
	if (NULL == xd)
		return exdrFILENOTFOUND;
	if ((result = do_trnheader(xd,1,&sh)) != exdrOK)
		return result;
	xdrfile_close(xd);
	*natoms = sh.natoms;
	
	return exdrOK;
}

int write_trr(XDRFILE *xd,int natoms,int step,float t,float lambda,
			  matrix box,rvec *x,rvec *v,rvec *f)
{
	return do_trn(xd,0,&step,&t,&lambda,box,&natoms,x,v,f);
}

int read_trr(XDRFILE *xd,int natoms,int *step,float *t,float *lambda,
			 matrix box,rvec *x,rvec *v,rvec *f)
{
	return do_trn(xd,1,step,t,lambda,box,&natoms,x,v,f);
}

}
}

