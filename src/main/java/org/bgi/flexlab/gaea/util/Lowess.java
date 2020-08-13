/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.util;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;

public class Lowess {
	private static class ShiftedLeft extends AbstractList<Double> {
		private List<Double> delegate;
		public ShiftedLeft(List<Double> delegate) {
			this.delegate=delegate;
		}
		
		@Override
		public Double set(int i, Double d) {
			if(i==0) throw new IllegalArgumentException();
			return this.delegate.set(i-1,d);
		}
		
		@Override
		public Double get(int i) {
			if(i==0) throw new IllegalArgumentException();
			return this.delegate.get(i-1);
		}

		@Override
		public int size() {
			return this.delegate.size()+1;
		}
		
	}
	

	/** proportion of points in the plot which influence the smooth at each value */
	private double smoother_span = 2.0/3.0;
	/** the number of iterations which should be performed. */
	private int nsteps = 3;
	/** used to speed up computation */
	private double delta_speed = -1;
	/** perform some basic checks */
	private boolean paranoid = true;
	
	private static class ListFactory<T> {
		public List<T> createList(int capacity)
			{
			return new ArrayList<T>(capacity);
			}
		public void disposeList(List<T> list)
			{
			if(list==null) return;
			list.clear();
			}
		public void close() {
			
			}
	}
	
	private  ListFactory<Double> listFactory = new ListFactory<Double>();
	
	public Lowess() {
	}
	
	public Lowess(double smoother_span) {
		this.smoother_span = smoother_span;
	}
	
	public Lowess(double smoother_span, int nsteps) {
		this.smoother_span = smoother_span;
		this.nsteps = nsteps;
	}
	
	public void setListFactory(ListFactory<Double> listFactory) {
		this.listFactory = listFactory;
	}
	
	private static double fsquare(double x) {
	    return x * x;
	}
	
	private static double fcube(double x) {
	    return x * x * x;
	}
	
	private static List<Double> shiftLeft(List<Double> L) {
		return new ShiftedLeft(L);
	}

	private static List<Double> subList(List<Double> L,int index) {
		return L.subList(index,L.size());
	}
	
	private static int TYPE_CMP(Double x,Double y,boolean nalast) {
		boolean nax = Double.isNaN(x), nay = Double.isNaN(y);
		    if (nax && nay)	return 0;
		    if (nax)		return nalast ? 1 : -1;
		    if (nay)		return nalast ? -1 : 1;
		    if (x < y)		return -1;
		    if (x > y)		return 1;
		    return 0;
	}
	
	private static void rPsort2(List<Double> x,int lo, int hi, int k) {
		    double v, w;
		///#define TYPE_CMP rcmp
	    boolean nalast=true;					
	    int L, R, i, j;					
									
	    for (L = lo, R = hi; L < R; ) {				
			v = x.get(k);						
			for(i = L, j = R; i <= j;) {				
			    while (TYPE_CMP(x.get(i), v, nalast) < 0) i++;		
			    while (TYPE_CMP(v, x.get(j), nalast) < 0) j--;		
			    if (i <= j) { w = x.get(i);x.set(i,x.get(j));i++;  x.set(j,w);--j; }
			}							
			if (j < k) L = i;					
			if (k < i) R = j;					
	    }
		//#undef TYPE_CMP
	}

	private static void rPsort(List<Double> L,int n,int k) {
		rPsort2(L, 0, n-1, k);
	}
	
	
	
	private void lowest(List<Double> x, List<Double> y, int n, List<Double> xs, List<Double> ys, int nleft, int nright,
		List<Double> w, boolean userw, List<Double> rw, boolean[] ok ) {
	    int nrt, j;
	    double a, b, c, h, h1, h9, r, range;

	    x=shiftLeft(x);
	    y=shiftLeft(y);
	    w=shiftLeft(w);
	    rw=shiftLeft(rw);

	    range = x.get(n)-x.get(1);
	    h = Math.max( xs.get(0)-x.get(nleft), x.get(nright)-xs.get(0) );
	    h9 = 0.999*h;
	    h1 = 0.001*h;

	    /* sum of weights */

	    a = 0.;
	    j = nleft;
	    while (j <= n) {

			/* compute weights */
			/* (pick up all ties on right) */
	
			w.set(j, 0.0);
			r = Math.abs(x.get(j) - xs.get(0));
			if (r <= h9) {
			    if (r <= h1)
			    	w.set(j,1.0);
			    else
			    	w.set(j, fcube(1.-fcube(r/h)));
			    if (userw)
			    		w.set(j, w.get(j) * rw.get(j));
			    a += w.get(j);
			}
			else if (x.get(j) > xs.get(0))
			    break;
			j = j+1;
	    }

	    /* rightmost pt (may be greater */
	    /* than nright because of ties) */

	    nrt = j-1;
	    if (a <= 0.)
	    	ok[0] = false;
	    else {
	    	ok[0] = true;

		/* weighted least squares */
		/* make sum of w[j] == 1 */

		for(j=nleft ; j<=nrt ; j++)
			w.set(j, w.get(j)/a);
		if (h > 0.) {
		    a = 0.;

		    /*  use linear fit */
		    /* weighted center of x values */

		    for(j=nleft ; j<=nrt ; j++)
		    	a += w.get(j) * x.get(j);
		    b = xs.get(0) - a;
		    c = 0.;
		    for(j=nleft ; j<=nrt ; j++)
		    	c += w.get(j)*fsquare(x.get(j)-a);
		    if (Math.sqrt(c) > 0.001*range) {
		    	b /= c;

			/* points are spread out */
			/* enough to compute slope */

			for(j=nleft; j <= nrt; j++)
			    w.set(j,w.get(j) * (b*(x.get(j)-a) + 1.0));
		    }
		}
		ys.set(0,0.0);
		for(j=nleft; j <= nrt; j++)
		    ys.set(0,ys.get(0) +  w.get(j) * y.get(j));
	    }
	}

	private void clowess( List<Double> x, List<Double> y, int n, double f, int nsteps, double delta,
		     List<Double> ys, List<Double> rw, List<Double> res) {
	    int i, iter, j, last, m1, m2, nleft, nright, ns;
	    boolean ok[]={true};
	    double alpha, c1, c9, cmad, cut, d1, d2, denom, r, sc;

	    if (n < 2) {
			ys.set(0,  y.get(0));
			return;
		}

	    /* nleft, nright, last, etc. must all be shifted to get rid of these: */
	    x=shiftLeft(x);
	    y=shiftLeft(y);
	    ys=shiftLeft(ys);

	    /* at least two, at most n points */
	    ns = Math.max(2, Math.min(n, (int)(f*n + 1e-7)));

	    /* robustness iterations */

	    iter = 1;
	    while (iter <= nsteps+1) {
			nleft = 1;
			nright = ns;
			last = 0;	/* index of prev estimated point */
			i = 1;		/* index of current point */
	
			for(;;) {
			    if (nright < n) {
					/* move nleft,  nright to right */
					/* if radius decreases */
		
					d1 = x.get(i) - x.get(nleft);
					d2 = x.get(nright+1) - x.get(i);
		
					/* if d1 <= d2 with */
					/* x[nright+1] == x[nright], */
					/* lowest fixes */
		
					if (d1 > d2) {
		
					    /* radius will not */
					    /* decrease by */
					    /* move right */
		
					    nleft++;
					    nright++;
					    continue;
					}
			    }
	
			    /* fitted value at x[i] */
			    lowest(subList(x,1), subList(y,1), n, subList(x,i), subList(ys,i), nleft, nright, res, iter>1, rw, ok);
			    if (!ok[0]) ys.set(i, y.get(i));
			    
			    /* all weights zero */
			    /* copy over value (all rw==0) */
			    if (last < i-1) {
					denom = x.get(i)-x.get(last);
		
					/* skipped points -- interpolate */
					/* non-zero - proof? */
		
					for(j = last+1; j < i; j++) {
					    alpha = (x.get(j)-x.get(last))/denom;
					    ys.set(j, alpha*ys.get(i) + (1.-alpha)*ys.get(last) );
					}
			    }
	
			    /* last point actually estimated */
			    last = i;
	
			    /* x coord of close points */
			    cut = x.get(last)+delta;
			    for (i = last+1; i <= n; i++) {
					if (x.get(i) > cut)
					    break;
					if (x.get(i) == x.get(last)) {
					    ys.set(i, ys.get(last));
					    last = i;
					}
			    }
			    i = Math.max(last+1, i-1);
			    if (last >= n)
			    	break;
			}
			/* residuals */
			for(i = 0; i < n; i++)
			    res.set(i, y.get(i+1) - ys.get(i+1));
	
			/* overall scale estimate */
			sc = 0.;
			for(i = 0; i < n; i++) sc += Math.abs(res.get(i));
			sc /= n;
	
			/* compute robustness weights */
			/* except last time */
	
			if (iter > nsteps)
			    break;
			/* Note: The following code, biweight_{6 MAD|Ri|}
			   is also used in stl(), loess and several other places.
			   --> should provide API here (MM) */
			for(i = 0 ; i < n ; i++)
			    rw.set(i,Math.abs(res.get(i)));
	
			/* Compute   cmad := 6 * median(rw[], n)  ---- */
			/* FIXME: We need C API in R for Median ! */
			m1 = n/2;
			/* partial sort, for m1 & m2 */
			rPsort(rw, n, m1);
			if(n % 2 == 0) {
			    m2 = n-m1-1;
			    rPsort(rw, n, m2);
			    cmad = 3.*(rw.get(m1)+rw.get(m2));
			}
			else { /* n odd */
			    cmad = 6.*rw.get(m1);
			}
		
			if(cmad < 1e-7 * sc) /* effectively zero */
			    break;
			c9 = 0.999*cmad;
			c1 = 0.001*cmad;
			for(i = 0 ; i < n ; i++) {
			    r = Math.abs(res.get(i));
			    if (r <= c1)
				rw.set(i,1.);
			    else if (r <= c9)
				rw.set(i,fsquare(1.-fsquare(r/cmad)));
			    else
				rw.set(i,0.);
			}
			iter++;
	    }
	}
	//see also ..R-2.11.0/src/library/stats/R/lowess.R

	public List<Double> lowess(List<Double> x, List<Double> y, int n) {
	    if(n<1) throw new IllegalArgumentException("n =" + n + "<1");
	    if(paranoid) {
			for(int i=0;i< n;++i) {
			    if(Double.isNaN(x.get(i))) throw new RuntimeException("NAN: x["+i+"]");
			    if(Double.isNaN(y.get(i))) throw new RuntimeException("NAN: y["+i+"]");
			    if(i>0){
					if(x.get(i-1)> x.get(i)) throw new RuntimeException("Data not sorted on i:" + i + "\t" + x.get(i-1) +" > " + x.get(i));
				}
			}
	    }
	    double delta=this.delta_speed;
	    if(delta<0.0) {
			delta=.01*(x.get(n-1)-x.get(0));
		}
	    List<Double> rw=this.listFactory.createList(n);
	    List<Double> ys=this.listFactory.createList(n);
	    List<Double> res=this.listFactory.createList(n);
	    for(int i=0;i< n;++i) {rw.add(0.0);ys.add(0.);res.add(0.0);}
	    
	    clowess(x, y, n, this.smoother_span, this.nsteps, delta, ys, rw, res);

	    this.listFactory.disposeList(rw);
	    this.listFactory.disposeList(res);
	    return ys;	
	}
	
}
