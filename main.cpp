/* 
 * File:   main.cpp
 * Author: Locyber <Locyber@cuttingplane.com>
 *
 * Created on 2015年9月21日, 上午9:21
 */

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>      // std::inner_product
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

#define OK_DONE         0
#define ERR_INVALIDARG  -1
#define ERR_NODIRECTOR  -2  
#define ERR_INVALIDDATA_MN  -3

typedef long long int64; 
typedef unsigned long long uint64;
typedef std::vector< std::vector<double> > vector2;

using namespace std;

// <editor-fold desc="模型相关常变量">
string  networkFile  = "GAResult_Network.txt";
string  resultFile  = "GAResult_Cost.txt";
string  workingPath = "/wwwroot/demo.ccpcsoft.com/data";
int DC,RC;
double  za      = 1.96;            // 标准正态离差
double  t       = 0.95;             // weight factor associated with the operating cost 运营成本系数
double  r       = 1.2;              // 零售商需求方差和均值的比值 即
double  beta    = 0.0144;        // weight factor associated with the transportation cost运输成本系数
double  theta   = 1.45;         // weight factor associated with the inventory cost库存成本系数
double  w       = 0.556;

vector<double> f;               // DC的固定选址成本
vector<double> a;               // DC的单位运输成本（从工厂到DC）
vector<double> L;               // DC的提前时间
vector<double> h;               // DC的单位库存持有成本
vector<double> k;               // DC的单位固定订货成本
vector<double> g;               // DC的固定运输成本（工厂到DC）

vector<double> mu;              // R需求
vector2 am;                     // alpha*mu
vector<double> c;               // coefficient of gamma(S)+H(S)

vector2 d;                      // 二维数组，距离
vector<double>  CJA;            //存储 C_{j,A(j)} 的值
// </editor-fold>

// <editor-fold desc="GA相关变量">
vector<int>     U;               // U 未分配零售商
vector<int>     O;               // O 已开仓库
vector< vector< int > > A;                       // A[j] 表示仓库j所服务的零售商 O[j]=1时A[j]里才能有存储     
vector<double>  rho;             //相当于每个零售商的平均费用
int     IT=0;                    // iteration
struct betalpha{
    double  value;
    int     order;
};
// </editor-fold>

// <editor-fold desc="时间相关常变量">
uint64 startTime;           // The time when GA execution starts
uint64 endTime;             // The end time of GA
uint64 rgaTime;             // The end time of RGA
// </editor-fold>

/**
 * 载入模型数据
 */
void loadModuleData(int DC,int RC,string fPath)
{   
    string folderPath = fPath;
    folderPath = folderPath.append("/var_t.txt");
    ifstream testSamplevt(folderPath.c_str());
    testSamplevt>>t;
    testSamplevt.close();
    
    folderPath = fPath;
    ifstream testSamplevr(folderPath.append("/var_r.txt").c_str());
    testSamplevr>>r;
    testSamplevr.close();
    
    folderPath = fPath;
    ifstream testSamplevbeta(folderPath.append("/var_beta.txt").c_str());
    testSamplevbeta>>beta;
    testSamplevbeta.close();
    
    folderPath = fPath;
    ifstream testSamplevtheta(folderPath.append("/var_theta.txt").c_str());
    testSamplevtheta>>theta;
    testSamplevtheta.close();
    
    folderPath = fPath;
    ifstream testSamplevw(folderPath.append("/var_w.txt").c_str());
    testSamplevw>>w;
    testSamplevw.close();
    
    folderPath = fPath;
    ifstream testSampleg(folderPath.append("/g.txt").c_str());
    double gvalue;
    for ( int i = 0; i < DC && testSampleg >> gvalue; i++) { g.push_back(gvalue);}
    testSampleg.close();
 
    folderPath = fPath;
    ifstream testSamplek(folderPath.append("/k.txt").c_str());
    double kvalue;
    for ( int i = 0; i < DC && testSamplek >> kvalue; i++) {k.push_back(kvalue);}
    testSamplek.close();
    
    folderPath = fPath;
    ifstream testSamplea(folderPath.append("/a.txt").c_str());
    double avalue;
    for ( int i = 0; i < DC && testSamplea >> avalue; i++) {a.push_back(avalue);}
    testSamplea.close();
    
    folderPath = fPath;
    ifstream testSampleL(folderPath.append("/L.txt").c_str());
    double Lvalue;
    for ( int i = 0; i < DC && testSampleL >> Lvalue; i++) {L.push_back(Lvalue);}
    testSampleL.close();
    
    folderPath = fPath;
    ifstream testSampleh(folderPath.append("/h.txt").c_str());
    double hvalue;
    for ( int i = 0; i < DC && testSampleh >> hvalue; i++) {h.push_back(hvalue);}
    testSampleh.close();
    
    folderPath = fPath;
    ifstream testSamplef(folderPath.append("/f.txt").c_str());
    double fvalue;
    for ( int i = 0; i < DC && testSamplef >> fvalue; i++) {f.push_back(fvalue);}
    testSamplef.close();
/*
    folderPath = fPath;
    ifstream testSamplec(folderPath.append("/c.txt").c_str());
    double cvalue;
    for ( int i = 0; i < DC && testSamplec >> cvalue; i++) {c.push_back(cvalue);}
    testSamplec.close();
 */
    for(int j=0;j<DC;j++)
    {
        c.push_back(sqrt(2*theta*h[j]*(k[j]+beta*g[j]))+theta*h[j]*za*sqrt(L[j]*r));
    }
    
    folderPath = fPath;
    ifstream testSamplemu(folderPath.append("/mu.txt").c_str());
    double muvalue;
    for ( int i = 0; i < RC && testSamplemu >> muvalue; i++) {mu.push_back(muvalue);}
    testSamplemu.close();
    
    folderPath = fPath;
    ifstream testSampled(folderPath.append("/d.txt").c_str());
    double dvalue;
    int row,col;
    vector<double> tmp;
    for ( int i = 0; i <RC*DC && testSampled >> dvalue; i++) 
    {
        row = i/DC;
        col = i%DC;        
        
        //cout<<i<<" :"<<row<<" "<<col<<" "<<dvalue<<endl;
        
        if(col == 0)   // 第一个 
        {
            tmp.clear();            
        } 
        tmp.push_back(dvalue);
        
        if(col==DC-1)
        {
            d.push_back(tmp);
        }
    
    }
    testSampled.close();
    
    
    for ( int i = 0; i <DC; i++) 
    {
        vector<double> amtmp;
        for(int j=0;j<RC;j++)
        {
            amtmp.push_back((a[j]+d[j][i])*mu[j]);
        }
        am.push_back(amtmp);
    }
    /*
    folderPath = fPath;
    ifstream testSampleam(folderPath.append("/am.txt"));
    double amvalue;
    row=0;
    col=0;
    vector<double> amtmp;
    for ( int i = 0; i <RC*DC && testSampleam >> amvalue; i++) 
    {
        row = i/RC;
        col = i%RC;        
        
        if(col == 0) 
        {
            amtmp.clear();            
        } 
        amtmp.push_back(amvalue);
        
        if(col == RC-1)
        {
            am.push_back(amtmp);
        }
    
    }
    testSampleam.close();  
    */
}

/**
 * 工具函数，计算数组中所有元素的和
 * @param a             输入的整形数组
 * @param num_elements  需要计算的长度
 * @return              整形，数组中指定长度的值的和
 */
int util_sumArray(vector< int > a, int num_elements)
{
   int i, sum=0;
   if(num_elements>a.size()) return sum;
   for (i=0; i<num_elements; i++)
   {
	 sum = sum + a[i];
   }
   return(sum);
}

int util_checkDirExists(string path)
{
    struct stat s;
    int err = stat(path.c_str(), &s);
    if(-1 == err) 
    {
        return -1;
    } else 
    {
        if(S_ISDIR(s.st_mode)) 
        {
            return 1;
        } 
        else 
        {
            return 0;
        }
    }
}

/**
 * 工具函数，获得数组中最小的元素
 * @param a             输入的数组       
 * @param count         需要计算的长度
 * @return              双精，数组中指定长度的值中最小的
 */
double util_getMinFromArr(vector< double> a,int count)
{
    double ret = a[0];
    if(count>a.size()) return 0;
    for(int i=1;i<count;i++)
    {
        ret = (ret<a[i])?ret:a[i];
    }
    return ret;
}

/**
 *  Returns the amount of milliseconds elapsed since the UNIX epoch. 
 * Works on both windows and linux. 
 */
uint64 GetTimeMs64()
{
    #ifdef _WIN32
     /* Windows */
     FILETIME ft;
     LARGE_INTEGER li;

     /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
      * to a LARGE_INTEGER structure. 
      */
     GetSystemTimeAsFileTime(&ft);
     li.LowPart = ft.dwLowDateTime;
     li.HighPart = ft.dwHighDateTime;

     uint64 ret = li.QuadPart;
     ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
     ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

     return ret;
    #else
     /* Linux */
     struct timeval tv;

     gettimeofday(&tv, NULL);

     uint64 ret = tv.tv_usec;
     /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
     ret /= 1000;

     /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
     ret += (tv.tv_sec * 1000);

     return ret;
    #endif
}

/**
 *  结构体betalpha的比较函数  
 */
bool structCompare(betalpha a, betalpha b)
{
    return a.value<b.value;
}

inline void getDCRC(int* DC, int* RC, string folderPath)
{   
    string tmp = folderPath;
    ifstream testSamplevt(tmp.append( "/var_m.txt").c_str());
    testSamplevt>>*DC;
    testSamplevt.close();
    
    tmp = folderPath;
    ifstream testSamplevr(tmp.append("/var_n.txt").c_str());
    testSamplevr>>*RC;
    testSamplevr.close();
}

void outputModel()
{
    ofstream testSamplevt("/home/outputModel/var_t.txt");
    testSamplevt.precision(16);
    testSamplevt<<t;
    testSamplevt.close(); 

    ofstream testSamplevr("/home/outputModel/var_r.txt");
    testSamplevr.precision(16);
    testSamplevr<<r;
    testSamplevr.close();  

    ofstream testSamplevbeta("/home/outputModel/var_beta.txt");
    testSamplevbeta.precision(16);
    testSamplevbeta<<beta;
    testSamplevbeta.close(); 

    ofstream testSamplevtheta("/home/outputModel/var_theta.txt");
    testSamplevtheta.precision(16);
    testSamplevtheta<<theta;
    testSamplevtheta.close(); 

    ofstream testSamplevw("/home/outputModel/var_w.txt");
    testSamplevw.precision(16);
    testSamplevw<<w;
    testSamplevw.close();

    ofstream testSamplef("/home/outputModel/f.txt");
    testSamplef.precision(16);
    for(int i = 0; i < DC; i++)
    {
        testSamplef << f[i] << endl;
    }
    testSamplef.close();

    ofstream testSamplec("/home/outputModel/c.txt");
    testSamplec.precision(16);
    for(int i = 0; i < DC; i++)
    {
         testSamplec << c[i] << endl;
    }
    testSamplec.close();

    ofstream testSamplemu("/home/outputModel/mu.txt");
    testSamplemu.precision(16);
    for(int i = 0; i < RC; i++)
    {
         testSamplemu << mu[i] << endl;
    }
    testSamplemu.close();

    ofstream testSampled("/home/outputModel/d.txt");
    testSampled.precision(16);
    for(int i = 0; i < RC; i++)
    {
        for(int j=0;j<DC;j++)
         {
             testSampled << d[i][j] << " ";
         }
         testSampled << endl;
     }
     testSampled.close();

    ofstream testSampleam("/home/outputModel/am.txt");
    testSampleam.precision(16);
    for(int i = 0; i < DC; i++)
    {
        for(int j=0;j<RC;j++)
         {
             testSampleam << am[i][j] << " ";
         }
         testSampleam << endl;
     }
    testSampleam.close();
}

int main(int argc,char* argv[]) 
{	   
    
    if(argc!=2)
    {
        return ERR_INVALIDARG;
    }
    string ID = argv[1];
    // todo: check ID   
    if (ID.length()!=32) return ERR_INVALIDARG;
    // todo: check ID path exits   
    string folderPath = workingPath.append("/");
    folderPath.append(ID);
    if(util_checkDirExists(folderPath)!=1) return ERR_NODIRECTOR;
      
    getDCRC(&DC,&RC,folderPath);
    if((DC<=0) || (RC<=0)) return ERR_INVALIDDATA_MN;
    
    loadModuleData(DC,RC,folderPath);
    // todo: checkModuleData
    
    
    int i,j,ii;
//    string resultPath = folderPath;
//    resultPath.append("/");
//    resultPath.append(resultFile);
//    ofstream cout;
//    cout.open(resultPath.c_str(),fstream::out);
//    cout.precision(8);
    
    outputModel();

    // Greedy Algorithm Starts
//    cout<<"=====  greedy algorithm  ====="<<endl;
//    startTime = GetTimeMs64();
//    cout<<"start time:"<<startTime<<endl;
  
    vector< vector< betalpha > > DCxRC;
    for(i=0;i<DC;i++)
    {
        vector<betalpha> btvector;
        for(j=0;j<RC;j++)
        {
            betalpha bt = {mu[j]*(d[j][i]+a[i]),j};
            btvector.push_back(bt);
        }
        sort(btvector.begin(),btvector.end(),structCompare);
        DCxRC.push_back(btvector);       
    }
  
    for(i=0;i<RC;i++)
    {
        U.push_back(1);
        rho.push_back(0);
    }
    for(j=0;j<DC;j++)
    {
        O.push_back(0);
        CJA.push_back(0);
        vector< int > Arow;
        for(i=0;i<RC;i++)
        {
            Arow.push_back(0);
        }
        A.push_back(Arow);
    }
   
    outputModel();

    while(true)
    {      
        int nn  = util_sumArray(U,RC);       
        if(nn==0) break;
        
        IT++;
        cout<<"iteration: "<<IT<<endl;
        
        double min=0;                   // 本次最优值
        vector< int > Sadd;             // 本次最优值对应的S*
        int J = 0;                      // 本次最优值对应的j	
        
        for(i=0;i<RC;i++)
        {
            Sadd.push_back(0);
        }
        
        for(j=0;j<DC;j++)
        {                            
            vector<int> ba;         
            for(i=0;i<RC;i++)       // for a selected DC, pick unchosen Retailers, and calc the transportation cost
            {                              
                if(U[DCxRC[j][i].order]==1)
                {
                    ba.push_back(DCxRC[j][i].order);
                }                
            }            

            for(i=0;i<nn;i++)   // nn = the total count of unchosen retailers
            {
                vector< int > S (RC,0);       
                vector< int > AS (RC,0);

                for(ii=0;ii<=i;ii++)    // 对于所有没有选中的retailer
                {
                    S[ba[ii]] = 1;    // 建一个S数组,把刚才排序好的靠前的ii个标记为1
                }
                for(ii=0;ii<RC;ii++)
                {
                    if(A[j][ii]+S[ii]>0)   
                    {
                        AS[ii]=1;           // 如果以前已经选择过，也标记1，                    
                    }
                }
               
                if(j+i==0)  //第一个计算得到的数值
                {                   
                    double x1 = (double)(std::inner_product(S.begin(),S.end(),am[j].begin(),0));
                    double x2 = (double)(std::inner_product(S.begin(),S.end(),mu.begin(),0));
                    min=(double)((f[j]+beta*x1+t*pow(x2,w)+c[j]*sqrt(x2))/x2);
                    
                    for(ii=0;ii<RC;ii++)
                    {
                        Sadd[ii]=S[ii];
                    }
                    J=j;
                }
                else
                {
                    double x01 = (double)(std::inner_product(AS.begin(),AS.end(),am[j].begin(),0));
                    double x02 = (double)(std::inner_product(AS.begin(),AS.end(),mu.begin(),0));
                    double x03 = (double)(std::inner_product(S.begin(),S.end(),mu.begin(),0));
                    double x04 = (double)((f[j]+beta*x01+t*pow(x02,w)+c[j]*sqrt(x02)-CJA[j])/x03);
                    if(x04<min)
                    {
                        min=x04;
                        for(ii=0;ii<RC;ii++)
                        {
                            Sadd[ii]=S[ii];
                        }
                        J=j;
                    }
                }
            }
        }
        O[J]=1;
        CJA[J]=min*std::inner_product(Sadd.begin(),Sadd.end(),mu.begin(),0)+CJA[J];
        for(i=0;i<RC;i++)
        {
            if(Sadd[i]==1)
            {
                rho[i]=min;
                U[i]=0;
                A[J][i]=1;
            }
        }

	cout<<"min:\t"<<min<<endl;
    }
  
//    endTime = GetTimeMs64();
//    cout<<"end time:"<<endTime<<endl;
//    cout<<"Open DCs number: "<<util_sumArray(O,DC)<<endl<<endl;
//    cout<<"时间（毫秒）: "<<startTime -endTime<<endl<<endl;

    string resultPath = folderPath;
    resultPath.append("/");
    resultPath.append(networkFile);
    ofstream cout;
    cout.open(resultPath.c_str(),fstream::out);
    cout.precision(8);
    
    for(j=0;j<DC;j++)
    {
        if(O[j]==1)
        {
            cout<<j<<" ";
            for(i=0;i<RC;i++)
            {
                if(A[j][i]==1)
                {
                    cout<<i<<" ";
                }
            }
            cout<<endl;
        }
    }
    cout.close();
    
    resultPath = folderPath;
    resultPath.append("/");
    resultPath.append(resultFile);
    ofstream fout;
    fout.open(resultPath.c_str(),fstream::out);
    fout.precision(8);
    
    fout<<"Iterations:"<<IT<<endl;
    double tcost = std::inner_product(rho.begin(),rho.end(),mu.begin(),0.0);
    fout<<"total cost: "<<tcost<<endl;


    vector< double > mp ;
    vector< double > op(DC,0);          // 存储每个warehouse对应的最优值
    //int FS[DC][RC] = {[0 ... DC-1][0 ... RC-1]=1};
    vector< vector < int > > FS;
    for(i=0;i<DC;i++)
    {
        vector < int > fsrow;
        for(j=0;j<RC;j++)
        {
            fsrow.push_back(1);
        }
        FS.push_back(fsrow);
    }
    
    for(i=0;i<RC;i++)
    {
        mp.push_back(mu[i]*rho[i]);
    }
    
    for(j=0;j<DC;j++)
    {
        double y1 = std::inner_product(FS[j].begin(),FS[j].end(),am[j].begin(),0);
        double y2 = std::inner_product(FS[j].begin(),FS[j].end(),mu.begin(),0);
        double y3 = std::inner_product(FS[j].begin(),FS[j].end(),mp.begin(),0);
        op[j] = (f[j]+beta*y1+t*pow(y2,w)+c[j]*sqrt(y2))/y3;
               
        //最优值由FS[j]决定
        for(;;)
        {
            int yes = 0;    //用于表明FS[j]是否发生改变,即theta是不是最优
            std::vector<betalpha> batArr;
            for(i=0;i<RC;i++)
            {
                double z1 = beta*(d[i][j]+a[j])-op[j]*rho[i];
                struct betalpha batV = {z1,i};
                batArr.push_back(batV);
            }
            vector< double > batm;
            for(i=0;i<RC;i++)
            {
                batm.push_back(batArr[i].value * mu[i]);
            }
            std::sort(batArr.begin(),batArr.end(),structCompare);
            

            for(i=0;i<RC;i++)
            {
                
                if(batArr[i].value<0)   //＜0的个数即是可能最优解的个数
                {                    
                    vector< int > zz (RC,0);
                    for(ii=0;ii<=i;ii++)
                    {
                        zz[batArr[ii].order]=1;
                    }
                                      
                    double rgax1 = std::inner_product(zz.begin(),zz.end(),batm.begin(),0);
                    double rgax2 = std::inner_product(zz.begin(),zz.end(),mu.begin(),0);
                    double rgax3 = std::inner_product(zz.begin(),zz.end(),mp.begin(),0);
                    double rgax4 = std::inner_product(zz.begin(),zz.end(),am[j].begin(),0);
                    double rgax5 = rgax1+t*pow(rgax2,w)+c[j]*sqrt(rgax2)+f[j];
                    double rgax6 = (f[j]+beta*rgax4+t*pow(rgax2,w)+c[j]*sqrt(rgax2))/rgax3;
                    
                    if(rgax5<0)
                    {
                        
                        if(rgax6<op[j])
                        {
                            op[j]=rgax6;
                            yes=1;
                            for(ii=0;ii<RC;ii++)
                            {
                                FS[j][ii]=zz[ii];
                            }
                        }	
                    }
                }
            }
            if(yes==0)
            {
                break;
            }
        }		   
    }
    rgaTime = GetTimeMs64();
    cout<<"rga time:"<<rgaTime<<endl;
    cout<<"rGA="<<1/util_getMinFromArr(op,DC)<<endl;
    cout<<"时间（毫秒）: "<<rgaTime-startTime<<endl;  
    cout.close();
      
    return OK_DONE;    
}
void leftmain()
{
/*    
    int i,j,ii;
    ofstream cout(workingPath);
    cout.precision(8);
    
    loadModuleData(DC,RC,);

    
    // Greedy Algorithm Starts
    cout<<"=====  greedy algorithm  ====="<<endl;
    startTime = GetTimeMs64();
    cout<<"start time:"<<startTime<<endl;
    
    //std::list <RCostArr> DCxRC(DC);
    //sCost DCxRC[DC][RC];
    betalpha DCxRC[DC][RC];
    for(i=0;i<DC;i++)
    {
        for(j=0;j<RC;j++)
        {
            //DCxRC[i][j] = {beta*(d[i][j]+a[j]),i};
            //DCxRC[i][j].value = beta*(d[j][i]+a[i]);
            //DCxRC[i][j].value = d[j][i]+a[i];
            DCxRC[i][j].value = mu[j]*(d[j][i]+a[i]);
            //DCxRC[i][j].dc = i;
            DCxRC[i][j].order = j;
        }
        sort(DCxRC[i],DCxRC[i]+RC,structCompare);
    }
    
    
    while(true)
    {      
        int nn  = util_sumArray(U,RC);       
        if(nn==0) break;
        
        IT++;
        //cout<<"iteration: "<<IT<<endl;
        
        double min=0;                   // 本次最优值
        int Sadd[RC] = {0};             // 本次最优值对应的S*
        int J = 0;                      // 本次最优值对应的j	
                 
        for(j=0;j<DC;j++)
        {                            
            vector<int> ba;         
            for(i=0;i<RC;i++)       // for a selected DC, pick unchosen Retailers, and calc the transportation cost
            {                              
                if(U[DCxRC[j][i].order]==1)
                {
                    ba.push_back(DCxRC[j][i].order);
                }                
            }            

            for(i=0;i<nn;i++)   // nn = the total count of unchosen retailers
            {
                int S[RC] = {0};          
                int AS[RC] = {0};

                for(ii=0;ii<=i;ii++)    // 对于所有没有选中的retailer
                {
                    S[ba[ii]] = 1;    // 建一个S数组,把刚才排序好的靠前的ii个标记为1
                    //S[DCxRC[j][ii].order] = 1;
                }
                for(ii=0;ii<RC;ii++)
                {
                    if(A[j][ii]+S[ii]>0)   
                    {
                        AS[ii]=1;           // 如果以前已经选择过，也标记1，                    
                    }
                }
               
                if(j+i==0)  //第一个计算得到的数值
                {                   
                    double x1 = (double)(std::inner_product(S,S+RC,am[j],0));
                    double x2 = (double)(std::inner_product(S,S+RC,mu,0));
                    min=(double)((f[j]+beta*x1+t*pow(x2,w)+c[j]*sqrt(x2))/x2);
                    
                    for(ii=0;ii<RC;ii++)
                    {
                        Sadd[ii]=S[ii];
                    }
                    J=j;
                }
                else
                {
                    double x01 = (double)(std::inner_product(AS,AS+RC,am[j],0));
                    double x02 = (double)(std::inner_product(AS,AS+RC,mu,0));
                    double x03 = (double)(std::inner_product(S,S+RC,mu,0));
                    double x04 = (double)((f[j]+beta*x01+t*pow(x02,w)+c[j]*sqrt(x02)-CJA[j])/x03);
                    if(x04<min)
                    {
                        min=x04;
                        for(ii=0;ii<RC;ii++)
                        {
                            Sadd[ii]=S[ii];
                        }
                        J=j;
                    }
                }
            }
        }
        O[J]=1;
        CJA[J]=min*std::inner_product(Sadd,Sadd+RC,mu,0)+CJA[J];
        for(i=0;i<RC;i++)
        {
            if(Sadd[i]==1)
            {
                rho[i]=min;
                U[i]=0;
                A[J][i]=1;
            }
        }
    }
    
    endTime = GetTimeMs64();
    cout<<"end time:"<<endTime<<endl;
    double tcost = std::inner_product(rho,rho+RC,mu,0.0);
    cout<<"total cost: "<<tcost<<endl;
    cout<<"Open DCs number: "<<util_sumArray(O,DC)<<endl<<endl;
    cout<<"时间（毫秒）: "<<startTime -endTime<<endl<<endl;

    for(j=0;j<DC;j++)
    {
        if(O[j]==1)
        {
            cout<<"WH "<<j<<": ";
            for(i=0;i<RC;i++)
            {
                if(A[j][i]==1)
                {
                    cout<<i<<", ";
                }
            }
            cout<<endl;
        }
    }

    cout<<"Iterations:"<<IT<<endl;


    double mp[RC];
    double op[DC];          // 存储每个warehouse对应的最优值
    int FS[DC][RC] = {[0 ... DC-1][0 ... RC-1]=1};

    for(i=0;i<RC;i++)
    {
        mp[i]=mu[i]*rho[i];
    }
    
    for(j=0;j<DC;j++)
    {
        double y1 = std::inner_product(FS[j],FS[j]+RC,am[j],0);
        double y2 = std::inner_product(FS[j],FS[j]+RC,mu,0);
        double y3 = std::inner_product(FS[j],FS[j]+RC,mp,0);
        op[j]=(f[j]+beta*y1+t*pow(y2,w)+c[j]*sqrt(y2))/y3;
               
        //最优值由FS[j]决定
        for(;;)
        {
            int yes = 0;    //用于表明FS[j]是否发生改变,即theta是不是最优
            std::vector<betalpha> batArr;
            for(i=0;i<RC;i++)
            {
                double z1 = beta*(d[i][j]+a[j])-op[j]*rho[i];
                struct betalpha batV = {z1,i};
                batArr.push_back(batV);
            }
            double batm[RC];
            for(i=0;i<RC;i++)
            {
                batm[i]=batArr[i].value * mu[i];
            }
            std::sort(batArr.begin(),batArr.end(),structCompare);
                                   
            for(i=0;i<RC;i++)
            {
                if(batArr[i].value<0)   //＜0的个数即是可能最优解的个数
                {                    
                    int zz[RC]={0};
                    for(ii=0;ii<=i;ii++)
                    {
                        zz[batArr[ii].order]=1;
                    }
                                      
                    double rgax1 = std::inner_product(zz,zz+RC,batm,0);
                    double rgax2 = std::inner_product(zz,zz+RC,mu,0);
                    double rgax3 = std::inner_product(zz,zz+RC,mp,0);
                    double rgax4 = std::inner_product(zz,zz+RC,am[j],0);
                    double rgax5 = rgax1+t*pow(rgax2,w)+c[j]*sqrt(rgax2)+f[j];
                    double rgax6 = (f[j]+beta*rgax4+t*pow(rgax2,w)+c[j]*sqrt(rgax2))/rgax3;
                    
                    if(rgax5<0)
                    {
                        if(rgax6<op[j])
                        {
                            op[j]=rgax6;
                            yes=1;
                            for(ii=0;ii<RC;ii++)
                            {
                                FS[j][ii]=zz[ii];
                            }	
                        }	
                    }
                }
            }
            if(yes==0)
            {
                break;
            }
        }		   
    }
    rgaTime = GetTimeMs64();
    cout<<"rga time:"<<rgaTime<<endl;
    cout<<"rGA="<<1/util_getMinFromArr(op,DC)<<endl;
    cout<<"时间（毫秒）: "<<rgaTime-startTime<<endl;
    
    cout.close();
    return 0;
 */ 
}


