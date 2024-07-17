#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iomanip> 
#include <algorithm> 
#include <vector> 
#include <map>
#include <functional>
#include <windows.h> 
#include "stat.h"
FILE* stream;

using namespace std;

//*******************************************************************
    template <class Function>
    __int64 time_call(Function&& f)  {
        __int64 begin = GetTickCount();
        f();
        return GetTickCount() - begin;
    }
    
      map<string, function<double(int*, int, int, int*)> >testcrit = {
      {"Kruskal",kruskalstatistic},{"Leman",lemanstatistic},{"Mood",moodstatistic},
      {"Fisher",fisherstatistic},{"Capon",caponstatistic},{"VanDerVarden",vandervardenstatistic},
      {"Klotz",klotzstatistic}
       };

 long int TestPerm(string critt, int* a, int kk, int n, int* m,vector<double>&h) {
        int k, j, l, r;
        double z;
        long int num;
        num=0;
        while(1 > 0) {
            j = n - 2;
        while (j >= 0 && a[j] >= a[j + 1]) j--;
        if (j < 0) return num;
        k = n - 1;
        while (a[j] >= a[k]) k--;
        swap1(a[j], a[k]);
        l = j + 1; r = n - 1;
        while (l < r)  swap1(a[l++], a[r--]);
        num++;
        z =testcrit[critt](a, kk, n, m);
        h.push_back(z);
      }
        return num;
    }

//#######################################################################

    void perm() {
        int n, * a, km, j,analitic;
        long int i, nnew;
        double s;

       
    vector<double>h;
    vector <double> wrange;
    vector <double> pw;
    long int nn,num;
    string crit,ss;
    int k, *astat;
    double *wrange1, *pw1;

    ifstream inf("omega.inp");
    ofstream out("omega.out");

    inf >> ss;
    inf >>crit;
    inf >> ss;
    inf >>k;
    inf >> ss;
    int* m = new int[k];
    for (int i = 0; i < k; i++)    inf >> m[i];
    inf.close();

        n = 0;
        for (i = 0; i < k; i++)  n += m[i];
        a = new int[n + 1];
        km = 0;
        for (i = 0; i < k; i++) {
            for (j = 0; j < m[i]; j++) a[j + km] = i + 1;
            km = km + m[i];
        }

        analitic=0;
        if (crit == "Wilcoxon") {
			analitic=1;
            nn = m[0] * m[1] + 1;
            wrange1 = new double[nn + 2];
            pw1 = new double[nn + 2];
            nn = wilcoxon_exact(m[0], m[1], wrange1, pw1);
        }
        if (crit == "Series") {
        	analitic=1;
            nn = m[0] * m[1] + 1;
            wrange1 = new double[nn + 2];
            pw1 = new double[nn + 2];
            nn = series_exact(m, wrange1, pw1);
        }
        if (crit == "Ansari") {
        	analitic=1;
            nn = 1 + int((m[0] * m[1]) / 2);
            astat = new int[nn];
            pw1 = new double[nn];
            wrange1 = new double[nn];
            nn = ansari_exact(m[0], m[1], astat, pw1);
            for (i = 0; i < nn; i++) wrange1[i] = double(astat[i]);
        }
			if(analitic==1){
				out<<crit << '\n';
     			out << nn << '\n';
    			for (int i = 0; i < k; i++)    out << m[i] << "   ";
     			out << '\n';
                 for (int i = 0; i < nn; i++)  out << (i + 1) << ":" << "  " << setprecision(12) << fixed << wrange1[i] << "      " << pw1[i] << '\n';
				out.close();
				return;
		}		
			  
        //Перестановка с повторением
      
       s =testcrit[crit](a, k, n, m);
        h.push_back(s);
        wcout << L" Permutation...";
        wcout << '\n';
        wcout << L" Samples: ";
        for (int i = 0; i <k; i++)    wcout <<m[i] << "   ";
        wcout << '\n';
        auto elapsed = time_call([&] {num=TestPerm(crit, a, k, n, m, h);});
        wcout << L" took "<< (elapsed / 1000.0) << L" s." << endl;
        wcout << '\n';
 
        wcout << L" Sorting...";
      
        num++;
        elapsed = time_call([&] {qsortRecursive(h.data(), num); });
        wcout << L" took " <<(elapsed/1000.0)<<L" s." << endl;

        h.push_back(-1);

        nn = 0; nnew = 0;
        for (i = 0; i < num; i++) {
            if (round(h[i] * 100000) / 100000 == round(h[static_cast<vector<double,allocator<double>>::size_type>(i) + 1] * 100000) / 100000) {
                nnew++;
            }
            else {
                pw.push_back(nnew+1);
                wrange.push_back(h[i]);
                nn += 1;
                nnew = 0;
            }
        }
        s = 0.;
        for (i = 0; i < nn; i++) {
            s = s + pw[i] / double(num);
            pw[i] = s;
        }
        
   
	 out<<crit << '\n';
     out << num << '\n';
     out << nn << '\n';
     for (int i = 0; i < k; i++)    out << m[i] << "   ";
     out << '\n';
     for (int i = 0; i < nn; i++)  out << (i + 1) << ":" << "  " << setprecision(12) << fixed << wrange[i] << "      " << pw[i] << '\n';
    out.close();
        
    }
//******************************************************************
int main() {
    perm();
    return 0;
}