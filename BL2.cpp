#include <string.h>
#include <bits/stdc++.h>
#include <math.h>
#define pb push_back
#define dc3 vector <tuple <double, double, double > >
#define d3 tuple <double, double, double>
#define mpiif map < pair <int, int>, double>
#define mpiiff map <pair <int, int> , pair < double, double > >
#define mpiivff map <pair <int, int>, vector <pair <double, double> > >
#define ff first
#define ss second
using namespace std;
double inf = (double)(1<<15);
double error = 0.0001;
double mx=0, mx2=0;

double MIN(double a, double b){
	return (a<b?a:b);
}

double MAX(double a, double b){
        return (a<b?b:a);
}

dc3 bl(vector <mpiif> &ups, vector < mpiif > &ds, vector <pair <int, int > > &argmin, double eps, double A, double B){	
	dc3 pos;
	double HMU, HMD, ZP, D, U;
	mpiivff hm;
	for (int i=0;i<(int)(ds.end()-ds.begin());i++){
		int a=0, b=0, zero=0;
		double z=inf;
		for (int jj=0;jj*eps<A && !zero;jj++){
			for (int kk=0;kk*eps<B && !zero;kk++){
				if (hm.find({jj, kk})!=hm.end()){
					for (pair <double, double> UD:hm[{jj, kk}]){
						if (UD.ss<z){
							int vl=1;
							for (auto tt=ups[i].begin();tt!=ups[i].end() && vl;tt++){
								if (jj+tt->ff.ff-argmin[i].ff<0 || ((double)(jj+tt->ff.ff-argmin[i].ff))*eps>=A || \
									kk+tt->ff.ss-argmin[i].ss<0 || ((double)(kk+tt->ff.ss-argmin[i].ss))*eps>=B) vl=0;
								else{
									if(hm.find({jj+tt->ff.ff-argmin[i].ff, kk+tt->ff.ss-argmin[i].ss})!=hm.end()){
										for (pair <double, double> UD2:hm[{jj+tt->ff.ff-argmin[i].ff, kk+tt->ff.ss-argmin[i].ss}]){
											HMU = UD2.ss;
											HMD = UD2.ff;
											ZP = UD.ss;
											D = ds[i][{tt->ff.ff, tt->ff.ss}];
											U = tt->ss;
											if((ZP+D<=HMD && ZP+U>=HMU) || (ZP+D>=HMD && ZP+D<HMU) || (ZP+U>HMD && ZP+U<=HMU)) vl=0;
										}
									}
								}
							}
							if (vl){
								a=jj;
								b=kk;
								z=UD.ss;
							}
						}
					}
				}else{
					int vl=1;
                                        for (auto tt=ups[i].begin();tt!=ups[i].end() && vl;tt++){
                                        	if (jj+tt->ff.ff-argmin[i].ff<0 || (jj+tt->ff.ff-argmin[i].ff)*eps>=A || \
                       			                        kk+tt->ff.ss-argmin[i].ss<0 || (kk+tt->ff.ss-argmin[i].ss)*eps>=B) vl=0;
                                                else{
                                                	if(hm.find({jj+tt->ff.ff-argmin[i].ff, kk+tt->ff.ss-argmin[i].ss})!=hm.end()){
								for (pair <double, double> UD2:hm[{jj+tt->ff.ff-argmin[i].ff, kk+tt->ff.ss-argmin[i].ss}]){
									HMU = UD2.ss;
                                                                	HMD = UD2.ff;
                                                                	D = ds[i][{tt->ff.ff, tt->ff.ss}];
                                                                	U = tt->ss;

                                        	        	        if((D<=HMD && U>=HMU) || (D>=HMD && D<HMU) || (U>HMD && U<=HMU)) vl=0;
                                                       	//upper height higher than minimum at region
							//down height lower than maximum at region-->middle
								}
							}
                                                }
                                        }
                                       	if (vl){
                         	                a=jj;
                                                b=kk;
                                                z=0.0;
						zero=1;
       	                                }
				}
			}
		}
		if (abs(z-inf)<=error){
			return pos;
		}else{

			for (auto tt= ups[i].begin();tt!=ups[i].end();tt++){
				double Zmn, Zmx;

				Zmn = z+ds[i][{tt->ff.ff, tt->ff.ss}];
				Zmx = z+tt->ss;
				
				if (hm.find({a+tt->ff.ff-argmin[i].ff, b+tt->ff.ss-argmin[i].ss})!=hm.end())hm[{a+tt->ff.ff-argmin[i].ff, b+tt->ff.ss-argmin[i].ss}].pb({Zmn, Zmx});
				else{
					vector <pair <double, double > > aux;
					aux.pb({Zmn, Zmx});
					hm[{a+tt->ff.ff-argmin[i].ff, b+tt->ff.ss-argmin[i].ss}] = aux;
				}
				mx = MAX(Zmx, mx);
			}
			pos.pb({(a-argmin[i].ff)*eps, (b-argmin[i].ss)*eps, z});
		}
		/*for (auto t=hm.begin();t!=hm.end();t++){
                                        cout << (t->ff.ff)*eps << ' ' << (t->ff.ss)*eps << endl;
                                        for (pair <double, double> e:t->ss){
                                                cout << '\t' << e.ff+mx2 << ' ' << e.ss+mx2 << endl;
                                        }
                                }
		cout << endl;*/
		//for (auto tt=hm.begin();tt!=hm.end();tt++)cout << tt->ff.ff << ' ' << tt->ff.ss << ' ' << tt->ss.ff << ' ' << tt->ss.ss << endl;
		//
	}
	/*for (auto t=hm.begin();t!=hm.end();t++){
                                        //cout << t->ff.ff << ' ' << t->ff.ss << endl;
                                        for (pair <double, double> e:t->ss){
                                                if (e.ff-e.ss>=0)cout << '\t' << e.ff << ' ' << e.ss << endl;
                                        }
                          	}*/

	return pos;
}

int main(){
	int nf;
	cin >> nf;
	vector < mpiif > ups, ds;
	vector <pair <int, int > > argmin;
	double eps, A, B;
	cin >> eps >> A >> B;
	string fname[nf];
	string hmname;
	for (int ii=0;ii<nf;ii++)cin >> fname[ii];
	cin >> hmname;
	
	int a, b, x, y;
	double z;
	for (int ii=0;ii<nf;ii++){
		std::ifstream ifs;
		ifs.open(fname[ii]);
		ifs >> a >> b;

		argmin.pb({a, b});
		ifs >> a;
		map <pair <int, int>, double> up;
		for (int i=0;i<a;i++){
			ifs >> x >> y >> z;
			up[{x, y}] = z;
		}
		ifs >> b;
		map <pair <int, int>, double> d;
		for (int i=0;i<b;i++){
			ifs >> x >> y >> z;
			d[{x, y}] = z;
		}
		/*for (auto t=up.begin();t!=up.end();t++){
			if (t->ss-d[{t->ff.ff, t->ff.ss}]<0){
				cout << fname[ii] << endl;
			}
		}*/
		ups.pb(up);
		ds.pb(d);
		ifs.close();
	}
	if (hmname.compare("default")!=0)1;
	dc3 pos;
	while (nf>0){
		mx=0;
		dc3 pos2 = bl(ups, ds, argmin, eps, A, B);
		int sz = (int)(pos2.end()-pos2.begin());
		//cout << sz << endl;
		nf-=sz;
		//cout << sz << endl;
		for (int i=0;i<sz;i++)pos.pb({get<0>(pos2[i]), get<1>(pos2[i]), get<2>(pos2[i])+mx2});
		for (int i=0;i<sz;i++)ups.erase(ups.begin());
		for (int i=0;i<sz;i++)ds.erase(ds.begin());
		for (int i=0;i<sz;i++)argmin.erase(argmin.begin());
		mx2+=mx;
	}
	
	int nobjs = (int)(pos.end()-pos.begin());
	cout << nobjs << endl;
	for (int i=0;i<nobjs;i++){
		cout << get<0>(pos[i]) << ' ' << get<1>(pos[i]) << ' ' << get<2>(pos[i]) << endl;
	}
        return 0;
}
