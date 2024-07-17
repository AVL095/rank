let s2pi = 2.50662827463100050242;

function invnormaldistribution(y0)  {

let expm2 = 0.13533528323661269189;
let maxrealnumber=Math.pow(10,300);
let minrealnumber=Math.pow(10,-300);

            let result = 0;
            let x = 0;
            let y = 0;
            let z = 0;
            let y2 = 0;
            let x0 = 0;
            let x1 = 0;
            let code = 0;
            let p0 = 0;
            let q0 = 0;
            let p1 = 0;
            let q1 = 0;
            let p2 = 0;
            let q2 = 0;
            let zz;
            
            if(y0<=0) {
                result = -maxrealnumber;
                return result;
            }
            if(y0>=1) {
                result = maxrealnumber;
                return result;
            }
            code = 1;
            y = y0;
            if(y>1.0-expm2) {
                y = 1.0-y;
                code = 0;
                }
            if(y>expm2) {
                y = y-0.5;
                y2 = y*y;
                p0 = -59.9633501014107895267;
                p0 = 98.0010754185999661536+y2*p0;
                p0 = -56.6762857469070293439+y2*p0;
                p0 = 13.9312609387279679503+y2*p0;
                p0 = -1.23916583867381258016+y2*p0;
                q0 = 1;
                q0 = 1.95448858338141759834+y2*q0;
                q0 = 4.67627912898881538453+y2*q0;
                q0 = 86.3602421390890590575+y2*q0;
                q0 = -225.462687854119370527+y2*q0;
                q0 = 200.260212380060660359+y2*q0;
                q0 = -82.0372256168333339912+y2*q0;
                q0 = 15.9056225126211695515+y2*q0;
                q0 = -1.18331621121330003142+y2*q0;
                x = y+y*y2*p0/q0;
                x = x*s2pi;
                result = x;
                return result;
            }
            zz=Math.log(y);
            x = Math.sqrt(-(2.0*zz));
            x0 = x-Math.log(x)/x;
            z = 1.0/x;
            if(x<8.0) {
                p1 = 4.05544892305962419923;
                p1 = 31.5251094599893866154+z*p1;
                p1 = 57.1628192246421288162+z*p1;
                p1 = 44.0805073893200834700+z*p1;
                p1 = 14.6849561928858024014+z*p1;
                p1 = 2.18663306850790267539+z*p1;
                p1 = -(1.40256079171354495875*0.1)+z*p1;
                p1 = -(3.50424626827848203418*0.01)+z*p1;
                p1 = -(8.57456785154685413611*0.0001)+z*p1;
                q1 = 1;
                q1 = 15.7799883256466749731+z*q1;
                q1 = 45.3907635128879210584+z*q1;
                q1 = 41.3172038254672030440+z*q1;
                q1 = 15.0425385692907503408+z*q1;
                q1 = 2.50464946208309415979+z*q1;
                q1 = -(1.42182922854787788574*0.1)+z*q1;
                q1 = -(3.80806407691578277194*0.01)+z*q1;
                q1 = -(9.33259480895457427372*0.0001)+z*q1;
                x1 = z*p1/q1;
            }
            else  {
                p2 = 3.23774891776946035970;
                p2 = 6.91522889068984211695+z*p2;
                p2 = 3.93881025292474443415+z*p2;
                p2 = 1.33303460815807542389+z*p2;
                p2 = 2.01485389549179081538*0.1+z*p2;
                p2 = 1.23716634817820021358*0.01+z*p2;
                p2 = 3.01581553508235416007*0.0001+z*p2;
                p2 = 2.65806974686737550832*0.000001+z*p2;
                p2 = 6.23974539184983293730*0.000000001+z*p2;
                q2 = 1;
                q2 = 6.02427039364742014255+z*q2;
                q2 = 3.67983563856160859403+z*q2;
                q2 = 1.37702099489081330271+z*q2;
                q2 = 2.16236993594496635890*0.1+z*q2;
                q2 = 1.34204006088543189037*0.01+z*q2;
                q2 = 3.28014464682127739104*0.0001+z*q2;
                q2 = 2.89247864745380683936*0.000001+z*q2;
                q2 = 6.79019408009981274425*0.000000001+z*q2;
                x1 = z*p2/q2;
            }
            x = x0-x1;
            if( code!=0 ) x = -x;
            result = x;
            return result;
          }
//********Нормальные порядковые статистики**************************************

function ordern(n,pr,ps) {

/*
  Computes David-Johnson approximation for mean and covariance between rth
  and sth order statistics from the normal dist. for a sample size n.
  pr=r/(n+1);ps=s/(n+1);
  vorder[0] - mean;
  vorder[1] - covariance between rth and sth order statistics 
  if r=s then vorder[1] - variance
*/
  let qr,qs,xr,xs,dr,ds;
  let xr1,xr2,xr3,xr4,xr5,xr6;
  let xs1,xs2,xs3,xs4,xs5,xs6;
  let z1,z2,z3,z4,z5,z6,z7;
  let er,crs;
  let vorder=[];

 qr=1-pr;qs=1-ps; 
xr=invnormaldistribution(pr);
xs=invnormaldistribution(ps);
dr=s2pi*Math.exp(xr*xr/2.);
ds=s2pi*Math.exp(xs*xs/2.);
xr1=dr;
xr2=xr*Math.pow(dr,2);
xr3=(2*xr*xr+1)*Math.pow(dr,3);
xr4=(6*xr*xr*xr+7*xr)*Math.pow(dr,4);
xr5=(24*Math.pow(xr,4)+46*xr*xr+7)*Math.pow(dr,5);
xr6=(120*Math.pow(xr,5)+326*Math.pow(xr,3)+127*xr)*Math.pow(dr,6);
xs1=ds;
xs2=xs*Math.pow(ds,2);
xs3=(2*xs*xs+1)*Math.pow(ds,3);
xs4=(6*Math.pow(xs,3)+7*xs)*Math.pow(ds,4);
xs5=(24*Math.pow(xs,4)+46*Math.pow(xs,2)+7)*Math.pow(ds,5);
xs6=(120*Math.pow(xs,5)+326*Math.pow(xs,3)+127*xs)*Math.pow(ds,6);

er=xr+pr*qr*xr2/(2*(n+2))+pr*qr*((qr-pr)*xr3/3+pr*qr*xr4/8)/Math.pow(n+2,2)+pr*qr*(-(qr-pr)*xr3/3+(Math.pow(qr-pr,2)-pr*qr)*xr4/4+qr*pr*(qr-pr)*xr5/6+Math.pow(qr*pr,2)*xr6/48)/Math.pow(n+2,3);

z1=(qr-pr)*xr2*xs1+(qs-ps)*xr1*xs2+pr*qr*xr3*xs1/2+ps*qs*xr1*xs3/2+pr*qs*xr2*xs2/2;
z1=z1*pr*qs/Math.pow(n+2,2);
z2=-(qr-pr)*xr2*xs1-(qs-ps)*xr1*xs2+(Math.pow(qr-pr,2)-pr*qr)*xr3*xs1;
z3=(Math.pow(qs-ps,2)-ps*qs)*xr1*xs3+(1.5*(qr-pr)*(qs-ps)+0.5*ps*qr-2*pr*qs)*xr2*xs2;
z4=(5/6)*pr*qr*(qr-pr)*xr4*xs1+(5/6)*ps*qs*(qs-ps)*xr1*xs4+(pr*qs*(qr-pr)+0.5*pr*qr*(qs-ps))*xr3*xs2;
z5=(pr*qs*(qs-ps)+0.5*ps*qs*(qr-pr))*xr2*xs3 + (1 /8)*Math.pow(pr*qr,2)*xr5*xs1+(1/8)*Math.pow(ps*qs,2)*xr1*xs5;
z6=0.25*Math.pow(pr,2)*qr*qs*xr4*xs2+0.25*pr*ps*Math.pow(qs,2)*xr2*xs4+(2*Math.pow(pr*qs,2)+3*pr*qr*ps*qs)*xr3*xs3/12;
z7=z2+z3+z4+z5+z6;
crs=z1+pr*qs*z7/Math.pow(n+2,3)+pr*qs*xr1*xs1/(n+2);

vorder[0]=er;
vorder[1]=crs;
return(vorder);
}
//***********cnm=n!/m!*(n-m)!********************************
function cnm(n,m) {
  var s1,s2,i;
  s1=0; s2=0;
  for (i=m+1;i<=n;i++) s1=s1+Math.log(i);
    for (i=1;i<=n-m;i++)   s2=s2+Math.log(i);
return Math.exp(s1-s2);
}



double normal_dist_inv_cdf(double p, double mu, double sigma) {

    double q,r,num,x,den;
  
    q = p - 0.5;
    if(fabs(q) <= 0.425) {
        r = 0.180625 - q * q;
        num = (((((((2.5090809287301226727e+3 * r +3.3430575583588128105e+4) * r + 6.7265770927008700853e+4) * r + 4.5921953931549871457e+4) * r +
                     1.3731693765509461125e+4) * r +1.9715909503065514427e+3) * r +1.3314166789178437745e+2) * r +3.3871328727963666080e+0) * q;
        den = (((((((5.2264952788528545610e+3 * r +2.8729085735721942674e+4) * r + 3.9307895800092710610e+4) * r + 2.1213794301586595867e+4) * r +
                     5.3941960214247511077e+3) * r + 6.8718700749205790830e+2) * r + 4.2313330701600911252e+1) * r +1.0);
        x = num / den;
        return(mu + (x * sigma));
     }

    if(q <= 0.0)  {
      r=p;
   }
    else {
      r=1.0-p;
    }
    r =sqrt(-log(r));
    if(r <= 5.0) {
        r = r - 1.6;
        num = (((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r +2.41780725177450611770e-1) * r +1.27045825245236838258e+0) * r +
                     3.64784832476320460504e+0) * r + 5.76949722146069140550e+0) * r +4.63033784615654529590e+0) * r +1.42343711074968357734e+0);
        den = (((((((1.05075007164441684324e-9 * r +5.47593808499534494600e-4) * r +1.51986665636164571966e-2) * r + 1.48103976427480074590e-1) * r +
                     6.89767334985100004550e-1) * r +1.67638483018380384940e+0) * r +2.05319162663775882187e+0) * r +
                     1.0);
    } else {
        r = r - 5.0;
        num = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r +1.24266094738807843860e-3) * r + 2.65321895265761230930e-2) * r +
                     2.96560571828504891230e-1) * r +1.78482653991729133580e+0) * r + 5.46378491116411436990e+0) * r + 6.65790464350110377720e+0);
        den = (((((((2.04426310338993978564e-15 * r +1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r +7.86869131145613259100e-4) * r +
                     1.48753612908506148525e-2) * r +1.36929880922735805310e-1) * r + 5.99832206555887937690e-1) * r +1.0);
    }
    x = num / den;
    if(q < 0.0)  x = -x;
    return(mu + (x * sigma));

}
