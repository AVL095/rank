#include <math.h>

double maxrealnumber=pow(10,100);
double minrealnumber=pow(10,-100);
double pi=3.14159265358979;
double spi=0.398942280401433;
double expm2 = 0.13533528323661269189;
double s2pi=2.50662827463100050242;
double invnormaldistribution(double y0);
void ordern(int n, double pr, double ps, double& er, double& vrs);
double cnm(int n, int m);
void qsortRecursive(double *mas, long int size);
double fisherstatistic(int* a,int kx,int n,int* m);
double vandervardenstatistic(int* a,int kx,int n,int* m);
double caponstatistic(int* a,int kx,int n,int* m);
double klotzstatistic(int* a,int kx,int n,int* m);
double moodstatistic(int* a, int kx,int n,int* m);
double lemanstatistic(int* a, int kx, int n, int* m);
double kruskalstatistic(int* aa, int kx, int n, int* m);
long int wilcoxon_exact(int m, int n, double *wrange, double *pw);
int series_exact(int* m, double *wrange,double *pw);
double ansaristatistic(int* a,int kx,int n, int* m);
double wilcoxonstatistic(int* a,int kx, int n,int* m);
void swap(double& ax, double& ay);
void swap1(int& ax, int& ay);
long int count_perm(int k, int* m);

int ansari_exact(int m, int n, int* astat, double* pw);
int* frqadd(int* f1, int* f2, int l1in, int l1out, int l2, int nstart);
int* imply(int* f1, int* f2, int l1in, int l1out, int l2, int noff);
int gscale(int test, int other, int* pw);
int start2(int n, int* f);
int start1(int n, int* f);
//#######################################################
double ansaristatistic(int* a,int kx,int n, int* m) {  //Ansari-Breadly
    int i, j, msmal;
    double r;
    msmal = m[0];
    r = 0;
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= msmal; j++)    if (a[i] == j) r += ((n*1.0 + 1.0)/2. - abs(i*1.0 - (n*1.0 + 1.0)/2.0));
    }
    return r;
}
//*********************************************
double wilcoxonstatistic(int* a,int kx, int n,int* m) {
    int msmal, i, j;
    double r;
    msmal = m[0];
    r = 0;
    for (i=1;i<= n;i++) {
        for (j=1;j<=msmal;j++) if (a[i]==j) r += i;
    }
    return r;
}
//################################################################
double invnormaldistribution(double y0) {
	double result, x, y, z, y2, x0, x1, p0, q0, p1, q1, p2, q2, zz;
	int code;

	if (y0 <= 0) {
		result = -maxrealnumber;
		return result;
	}
	if (y0 >= 1) {
		result = maxrealnumber;
		return result;
	}
	code = 1;
	y = y0;
	if (y > 1.0 - expm2) {
		y = 1.0 - y;
		code = 0;
	}
	if (y > expm2) {
		y = y - 0.5;
		y2 = y * y;
		p0 = -59.9633501014107895267;
		p0 = 98.0010754185999661536 + y2 * p0;
		p0 = -56.6762857469070293439 + y2 * p0;
		p0 = 13.9312609387279679503 + y2 * p0;
		p0 = -1.23916583867381258016 + y2 * p0;
		q0 = 1;
		q0 = 1.95448858338141759834 + y2 * q0;
		q0 = 4.67627912898881538453 + y2 * q0;
		q0 = 86.3602421390890590575 + y2 * q0;
		q0 = -225.462687854119370527 + y2 * q0;
		q0 = 200.260212380060660359 + y2 * q0;
		q0 = -82.0372256168333339912 + y2 * q0;
		q0 = 15.9056225126211695515 + y2 * q0;
		q0 = -1.18331621121330003142 + y2 * q0;
		x = y + y * y2 * p0 / q0;
		x = x * s2pi;
		result = x;
		return result;
	}

	zz = log(y);
	x = sqrt(-(2.0 * zz));
	x0 = x - log(x) / x;
	z = 1.0 / x;

	if (x < 8.0) {
		p1 = 4.05544892305962419923;
		p1 = 31.5251094599893866154 + z * p1;
		p1 = 57.1628192246421288162 + z * p1;
		p1 = 44.0805073893200834700 + z * p1;
		p1 = 14.6849561928858024014 + z * p1;
		p1 = 2.18663306850790267539 + z * p1;
		p1 = -(1.40256079171354495875 * 0.1) + z * p1;
		p1 = -(3.50424626827848203418 * 0.01) + z * p1;
		p1 = -(8.57456785154685413611 * 0.0001) + z * p1;
		q1 = 1;
		q1 = 15.7799883256466749731 + z * q1;
		q1 = 45.3907635128879210584 + z * q1;
		q1 = 41.3172038254672030440 + z * q1;
		q1 = 15.0425385692907503408 + z * q1;
		q1 = 2.50464946208309415979 + z * q1;
		q1 = -(1.42182922854787788574 * 0.1) + z * q1;
		q1 = -(3.80806407691578277194 * 0.01) + z * q1;
		q1 = -(9.33259480895457427372 * 0.0001) + z * q1;
		x1 = z * p1 / q1;
	}
	else {
		p2 = 3.23774891776946035970;
		p2 = 6.91522889068984211695 + z * p2;
		p2 = 3.93881025292474443415 + z * p2;
		p2 = 1.33303460815807542389 + z * p2;
		p2 = 2.01485389549179081538 * 0.1 + z * p2;
		p2 = 1.23716634817820021358 * 0.01 + z * p2;
		p2 = 3.01581553508235416007 * 0.0001 + z * p2;
		p2 = 2.65806974686737550832 * 0.000001 + z * p2;
		p2 = 6.23974539184983293730 * 0.000000001 + z * p2;
		q2 = 1;
		q2 = 6.02427039364742014255 + z * q2;
		q2 = 3.67983563856160859403 + z * q2;
		q2 = 1.37702099489081330271 + z * q2;
		q2 = 2.16236993594496635890 * 0.1 + z * q2;
		q2 = 1.34204006088543189037 * 0.01 + z * q2;
		q2 = 3.28014464682127739104 * 0.0001 + z * q2;
		q2 = 2.89247864745380683936 * 0.000001 + z * q2;
		q2 = 6.79019408009981274425 * 0.000000001 + z * q2;
		x1 = z * p2 / q2;
	}
	x = x0 - x1;
	if (code != 0) x = -x;
	result = x;
	return result;
}
//**********************************************************************
void ordern(int n, double pr, double ps, double& er, double& vrs) {
    double p, pr1, ps1, xr, xr1, xr2, xr3, xr4, xr5, xr6, dr, qr, qs;
    double xs1, xs2, xs3, xs4, xs5, xs6, ds, xs;
    double z1, z2, z3, z4, z5, z6, z7;

    p = 1;
    xr = invnormaldistribution(pr);
    xs = invnormaldistribution(ps);
    qr = 1. - pr;
    qs = 1. - ps;
    pr1 = pr * p; ps1 = ps * p;
    xr = invnormaldistribution(pr1);
    xs = invnormaldistribution(ps1);
    dr = spi * exp(-xr * xr / 2.);
    ds = spi * exp(-xs * xs / 2.);
    xr1 = p / dr; xr2 = xr * xr1 * xr1;
    xr3 = (2. * xr * xr + 1.) * pow((p / dr), 3);
    xr4 = (6. * xr * xr * xr + 7. * xr) * pow((p / dr), 4);
    xr5 = (24. * pow(xr, 4) + 46. * xr * xr + 7.) * pow((p / dr), 5);
    xr6 = (120. * pow(xr, 5) + 326. * xr * xr * xr + 127. * xr) * pow((p / dr), 6);
    xs1 = p / ds; xs2 = xs * xs1 * xs1;
    xs3 = (2. * xs * xs + 1.) * pow((p / ds), 3);
    xs4 = (6. * xs * xs * xs + 7. * xs) * pow((p / ds), 4);
    xs5 = (24. * pow(xs, 4) + 46. * xs * xs + 7.) * pow((p / ds), 5);
    xs6 = (120. * pow(xs, 5) + 326. * xs * xs * xs + 127. * xs) * pow(p / ds, 6);

    er = xr + pr * qr * xr2 / (2. * (n + 2.)) + pr * qr * ((qr - pr) * xr3 / 3. + pr * qr * xr4 / 8.) / pow((n + 2.), 2) + pr * qr * (-(qr - pr) * xr3 / 3. + (pow((qr - pr), 2) - pr * qr) * xr4 / 4. + qr * pr * (qr - pr) * xr5 / 6. + pow((qr * pr), 2) * xr6 / 48.) / pow((n + 2.), 3);

    z1 = (qr - pr) * xr2 * xs1 + (qs - ps) * xr1 * xs2 + pr * qr * xr3 * xs1 / 2. + ps * qs * xr1 * xs3 / 2. + pr * qs * xr2 * xs2 / 2.;
    z1 = z1 * pr * qs / pow((n + 2.), 2);
    z2 = -(qr - pr) * xr2 * xs1 - (qs - ps) * xr1 * xs2 + (pow((qr - pr), 2) - pr * qr) * xr3 * xs1;
    z3 = (pow((qs - ps), 2) - ps * qs) * xr1 * xs3 + (1.5 * (qr - pr) * (qs - ps) + 0.5 * ps * qr - 2. * pr * qs) * xr2 * xs2;
    z4 = (5. / 6.) * pr * qr * (qr - pr) * xr4 * xs1 + (5. / 6.) * ps * qs * (qs - ps) * xr1 * xs4 + (pr * qs * (qr - pr) + .5 * pr * qr * (qs - ps)) * xr3 * xs2;
    z5 = (pr * qs * (qs - ps) + 0.5 * ps * qs * (qr - pr)) * xr2 * xs3 + (1. / 8.) * pow((pr * qr), 2) * xr5 * xs1 + (1. / 8.) * pow((ps * qs), 2) * xr1 * xs5;
    z6 = 0.25 * pr * pr * qr * qs * xr4 * xs2 + 0.25 * pr * ps * qs * qs * xr2 * xs4 + (2. * (pr * pr * qs * qs) + 3. * pr * qr * ps * qs) * xr3 * xs3 / 12.;
    z7 = z2 + z3 + z4 + z5 + z6;
    vrs = z1 + pr * qs * z7 / pow((n + 2.), 3) + pr * qs * xr1 * xs1 / (n + 2.);
}
//***********cnm=n!/m!*(n-m)!********************************
double cnm(int n, int m) {
    double s1, s2;
    int i;
    s1 = 0; s2 = 0;
    for (i = m + 1; i <= n; i++) s1 += log(i);
    for (i = 1; i <= n - m; i++) s2 += log(i);
    return exp(s1 - s2);
}
//*************************************************************
void qsortRecursive(double *mas, long int size) {
    long int i = 0;
    long int j = size - 1;
    double tmp,mid;
    mid = mas[size / 2];

    do {
        while (mas[i] < mid)  i++;
        while (mas[j] > mid)  j--;
        if (i <= j) {
            tmp = mas[i];mas[i] = mas[j];mas[j] = tmp; i++; j--;
        }
    } while (i <= j);
    if (j > 0)  qsortRecursive(mas, j + 1);
    if (i < size) qsortRecursive(&mas[i], size - i);
}
//*************************************************************
    double fisherstatistic(int* a,int kx,int n,int* m) {
      int i;
      double er,r,vrs;
      r=0;
       for (i=0;i<n;i++) { 
         ordern(n,(i+1.)/(n+1.),(i+1.)/(n+1.),er,vrs);
         if (a[i]==1) r+=er;
        }
         return r; 
     }
//**************************************************************
    double vandervardenstatistic(int* a,int kx,int n,int* m) {
      int i;
      double er,r;
      r=0;
       for (i=0;i<n;i++) { 
         er=invnormaldistribution((i+1.)/(m[0]+m[1]+1.));
         if (a[i]==1) r+=er;
        }
         return r; 
     }
//**************************************************************
    double caponstatistic(int* a,int kx,int n,int* m) {
      int i;
      double er,r,vrs;
      r=0;
       for (i=0;i<n;i++) { 
        ordern(n,(i+1.)/(n+1.),(i+1.)/(n+1.),er,vrs);
        er=vrs+er*er;
        if (a[i]==1) r+=er;
       }
         return r; 
     }
//**************************************************************
 double klotzstatistic(int* a,int kx,int n,int* m) {
      int i;
      double er,r;
      r=0;
       for (i=0;i<n;i++) { 
          er=invnormaldistribution((i+1.)/(m[0]+m[1]+1.));
          if (a[i]==1) r+=er*er;
        }
        
   return r;
  }
//********************************************************
double moodstatistic(int* a, int kx,int n,int* m) {
    int i, msmal;
    double r;
    msmal = 2;
    if (m[0] < m[1]) msmal = 1;
    r= 0;
    for (i = 0; i <n; i++)  if (a[i]==msmal) r+=pow((i+1.- (n + 1.) / 2.),2);
    return r;
}
//******************************************************************
double lemanstatistic(int* a, int kx, int n, int* m) {
    int i, k1,k2;
    double zleman,r1,r2;
    r1 = 0.; r2 = 0.; k1 = 0; k2 = 0;
    for (i = 0; i < n; i++) {
        if (a[i] == 1) {
            r1=r1+ (i - k1) * (i - k1); k1++;
        }
        if (a[i] == 2) {
            r2 += (i - k2) * (i - k2); k2++;
        }
    }
    zleman = (r1 * m[0] + r2 * m[1] + m[0] * m[1] / 6.) / (m[0] * m[0] * m[1] * m[1]) - 2. / 3.;
    return zleman;
}
//**********************************************************************
double kruskalstatistic(int* aa, int kx, int n, int* m) {
    int i,j;
    double s,r;
    s=0.;
     for (j=0;j<kx;j++) {
        r=0.;
        for (i=0;i<n;i++) if (aa[i]==j+1) r+=i+1;
        s=s+r*r/m[j];
    }
    return 12.*s/(1.0*n*(n+1.))-3.*(1.0*n+1.0);
}
//*********************Ansari**start1********************************************
int start1(int n, int* f) {
    int i, lout;
    lout = int(1 + n / 2);
    for (i = 1; i <= lout; i++) f[i] = 2;
    if ((n % 2) == 0) f[lout] = 1;
    return(lout);
}
//**********************Ansari start2**************************************
int start2(int n, int* f) {
    int one, two, three, four, i, j, a, b, lt1, ndo, nu, lout;

    one = 1; two = 2; three = 3; four = 4;
    nu = n - n % 2;
    j = nu + 1; lout = j; lt1 = lout + 1;
    ndo = int(lt1 / 2);
    a = one; b = three;
    for (i = 1; i <= ndo; i++) {
        f[i] = a; f[j] = a; j = j - 1; a = a + b; b = four - b;
    }
    if (nu == n) return(lout);
    nu = ndo + 1;
    for (i = nu; i <= lout; i++) f[i] = f[i] + two;
    f[lt1] = two; lout = lt1;
    return(lout);
}
//***************Ansari frqadd***********************************************
int* frqadd(int* f1, int* f2, int l1in, int l1out, int l2, int nstart) {

    int i1, i2, nxt, *fadd;
    
     fadd = new int[3];
    i2 = 1;
    for (i1 = nstart; i1 <= l1in; i1++) {
        f1[i1] = f1[i1] + 2 * f2[i2];
        i2++;
    }
    nxt = l1in + 1;
    l1out = l2 + nstart - 1;
    for (i1 = nxt; i1 <= l1out; i1++) {
        f1[i1] = 2 * f2[i2];
        i2++;
    }
    nstart++;
    fadd[1] = l1out; fadd[2] = nstart;
    return(fadd);
}
//***************Ansari imply*************************************************
int* imply(int* f1, int* f2, int l1in, int l1out, int l2, int noff) {
    int sum, diff, i2, i1, j2, j1, j2min, ndo;
	int *fimply;

     fimply = new int[4];
     
    i2 = 1 - noff; j1 = l1out; j2 = l1out - noff; l2 = j2;
    j2min = int((j2 + 1) / 2);
    ndo = int((l1out + 1) / 2);

    for (i1 = 1; i1 <= ndo; i1++) {
        if (i2 > 0) {
            sum = f1[i1] + f2[i2];
            f1[i1] = sum;
        }
        else {
            sum = f1[i1];
        }
        i2 = i2 + 1;
        if (j2 >= j2min) {
            if (j1 <= l1in) {
                diff = sum - f1[j1];
            }
            else {
                diff = sum;
            }

            f2[i1] = diff; f2[j2] = diff; j2 = j2 - 1;
        }
        f1[j1] = sum; j1 = j1 - 1;
    }
    fimply[1] = l1out; fimply[2] = l2; fimply[3] = noff;
    return(fimply);
}
//**************Ansari gscale*******************************************
int gscale(int test, int other, int* pw) {

    int  i, m, lres, mm1, nm1, nm2, ier, mnow, ks, j, ndo;
    int n, ln1, ln2, nc, l1out, l2out, n2b1, n2b2, ln3, kk, ai, z;
    int* a2, * a3, * fadd, * fimply;
    bool symm;

    ln1 = 0; ln2 = 0; l1out = 0; l2out = 0; ln1 = 0; ln2 = 0; ln3 = 0; n2b1 = 0; n2b2 = 0; kk = 0; nc = 0;
    ndo = 0; ier = 0; ks = 0; j = 0; i = 0;

    m = int(fmin(test, other));
    if (m < 0) return(0);
    n = int(fmax(test, other));
    lres = 1 + int((m * n) / 2);
    fadd = new int[lres];
    fimply = new int[lres];
    a2 = new int[lres];
    a3 = new int[lres];

    symm = false;
    z = (m + n) % 2;
    if (z == 0) symm = true;
    mm1 = m - 1;
    //*****************************************************         
    if (m <= 2) {
        if (mm1 < 0) {
            pw[1] = 1; return(lres);
        }

        if (mm1 == 0) ln1 = start1(n, pw);
        if (mm1 > 0)  ln1 = start2(n, pw);
        if (symm || (other > test)) return(lres);
        j = lres;
        ndo = int(lres / 2);
        for (i = 1; i <= ndo; i++) {
            ai = pw[i];
            pw[i] = pw[j];
            pw[j] = ai;
            j = j - 1;
        }
        return(lres);
    }
    //***********************************************************
    nm1 = n - 1; nm2 = n - 2; mnow = 3; nc = 3; ier = 0;
    //**************************************************************
    while (true) {
        if (ier == 0) {
            if ((n % 2) != 1) {
                n2b1 = 3;
                n2b2 = 2;
                ln1 = start2(n, pw);
                ln3 = start2(nm2, a3);
                ln2 = start1(nm1, a2);
                //***********************************************************************          
                fadd = frqadd(a2, a3, ln2, l2out, ln3, n2b2);
                l2out = fadd[1]; n2b2 = fadd[2];
                ln2 = ln2 + nm1;
                fimply = imply(a2, a3, l2out, ln2, j, nc);
                ln2 = fimply[1]; j = fimply[2]; nc = fimply[3];
                //***************************************************************************
                nc = nc + 1;
                if (mnow == m) break;
                mnow = mnow + 1;
            }
            else {
                n2b1 = 2;
                n2b2 = 3;
                ln1 = start1(n, pw);
                ln2 = start2(nm1, a2);
            }
        }
        //*******************************************************************
        fadd = frqadd(pw, a2, ln1, l1out, ln2, n2b1);
        l1out = fadd[1]; n2b1 = fadd[2];
        ln1 = ln1 + n;
        fimply = imply(pw, a3, l1out, ln1, ln3, nc);
        ln1 = fimply[1]; ln3 = fimply[2]; nc = fimply[3];
        nc = nc + 1;
        if (mnow == m) break;
        mnow = mnow + 1;
        fadd = frqadd(a2, a3, ln2, l2out, ln3, n2b2);
        l2out = fadd[1]; n2b2 = fadd[2];
        ln2 = ln2 + nm1;
        fimply = imply(a2, a3, l2out, ln2, j, nc);
        ln2 = fimply[1]; j = fimply[2]; nc = fimply[3];
        //***************************************************************************
        nc = nc + 1;
        if (mnow == m) break;
        mnow = mnow + 1;
        ier = 1;
    }
    //********************************************************************
    if (symm) return(lres);
    ks = int((m + 3) / 2);
    j = 1;
    for (i = ks; i <= lres; i++) {
        if (i > ln1) {
            pw[i] = a2[j];
        }
        else {
            pw[i] = pw[i] + a2[j];
        }
        j = j + 1;
    }
    if (other < test) return(lres);
    j = lres;
    ndo = int(lres / 2);
    for (i = 1; i <= ndo; i++) {
        ai = pw[i];
        pw[i] = pw[j];
        pw[j] = ai;
        j = j - 1;
    }
    return(lres);
}
//************************
int ansari_exact(int m, int n, int* astat, double* pw) {
    int min_val, max_val, i, nrows, a0;
    double sum, s;
    int* w;

    min_val = m;
    max_val = n;
    nrows = 1 + m * n;
    w = new int[nrows];
    
    nrows = gscale(min_val, max_val, w);
    a0 = int((min_val + 1) / 2) * (1 + int(min_val / 2));
    
     sum = 0;
    for (i = 1; i <= nrows; i++) {
        astat[i - 1] = a0 + i - 1;
        sum += w[i];
    }
    
     s = 0;
    for (i = 0; i < nrows; i++) {
        s = s + w[i+1];
        pw[i] = double(s/sum);
    }
    
    delete[] w;
    return(nrows);
}
//************ Р С™РЎР‚Р С‘РЎвЂљР ВµРЎР‚Р С‘Р в„– РЎРѓР ВµРЎР‚Р С‘Р в„– РЎвЂљР С•РЎвЂЎР Р…Р С•Р Вµ РЎР‚Р В°РЎРѓР С—РЎР‚Р ВµР Т‘Р ВµР В»Р ВµР Р…Р С‘Р Вµ**(Р В°Р Р…Р В°Р В»Р С‘РЎвЂљР С‘Р С”Р В°)************
int series_exact(int* m, double *wrange,double *pw) {
    int nn,k, kcur, m1, m2, n, i, ksr;
    double pc, w2, s;
    double *w;

    nn = m[0] * m[1] + 1;
    w = new double[nn];

    k = 2;
    m1 = m[0]; m2 = m[1]; n = 0;
    for (i = 0; i < k; i++)  n = n + m[i];
    w2 = cnm(n, m1);
    s = 0;
    for (i = 1; i <= 5000; i++) {
        kcur = i + 1;
        if (kcur % 2 == 0) {
            pc = double(2. * cnm(m1 - 1, kcur / 2 - 1) * cnm(m2 - 1, kcur / 2 - 1));
        }
        else {
            pc = double(cnm(m1 - 1, (kcur - 1) / 2) * cnm(m2 - 1, (kcur - 1) / 2 - 1) + cnm(m1 - 1, (kcur - 1) / 2 - 1) * cnm(m2 - 1, (kcur - 1) / 2));
        }
        s += pc;
        w[i-1]=pc;
        pw[i-1]=(s / w2);
        wrange[i-1]=double(kcur);
        if (pw[i - 1] > 1) break;
    }
    ksr = i - 1;
    return ksr;
}
//********Р С™РЎР‚Р С‘РЎвЂљР ВµРЎР‚Р С‘Р в„– Р Р€Р С‘Р В»Р С”Р С•Р С”РЎРѓР С•Р Р…Р В° (РЎвЂљР С•РЎвЂЎР Р…Р С•Р Вµ РЎР‚Р В°РЎРѓР С—РЎР‚Р ВµР Т‘Р ВµР В»Р ВµР Р…Р С‘Р Вµ)*(Р В°Р Р…Р В°Р В»Р С‘РЎвЂљР С‘Р С”Р В°)***********************

long int wilcoxon_exact(int m, int n, double *wrange, double *pw) {

    int* work, minmn, maxmn, inx;
    int* w;
    int n1, kk, k, j, i,sm;
    long int mn;
    double sum;

    mn = m * n + 1; maxmn = int(fmax(m, n)); minmn =int(fmin(m, n)); n1 = maxmn + 1;
    work = new int[mn + 2]; w = new int[mn + 2];
   
    for (i = 1; i <= n1; i++)  w[i] = 1;
    n1++;
    for (i = n1; i <= mn; i++) w[i] = 0;
    work[1] = 0; inx = maxmn;

    for (i = 2; i <= minmn; i++) {
        work[i] = 0; inx = inx + maxmn; n1 = inx + 2; kk = 1 + inx / 2; k = i;
        for (j = 1; j <= kk; j++) {
            k++;  n1--; sm = w[j] + work[j];
            w[j] = sm; work[k] = sm - w[n1]; w[n1] = sm;
        }
    }

    for (i = 0; i < mn; i++)  pw[i]=double(w[i + 1]);
    sum = 0;
    for (i = 0; i < mn; i++) {
        sum += pw[i];
        wrange[i]=double(i + minmn * (minmn + 1.) / 2.); //Wilcoxon
        //wrange.push_back((double(i)); //Mann-Whitney
        pw[i] = sum;
    }
    for (i = 0; i < mn; i++) pw[i] = pw[i] / sum;
    delete[] work;
    return mn;
}
//************************************************************
void swap1(int& ax, int& ay) {
        int s = ax; ax = ay; ay = s;
    }
void swap(double& ax, double& ay) {
            double s = ax; ax = ay; ay = s;
        }
//********************************************************
long int count_perm(int k, int* m) {
        int i, n;
        long int knum;
        double z, w;
        w = 1.; n = 0;
        for (i = 0; i < k; i++) {
            z = w * tgamma((m[i]) + 1.); w = z; n = n + m[i];
        }
        w = tgamma(n + 1.) / w;
        knum = long(w);
        return knum;
    }
