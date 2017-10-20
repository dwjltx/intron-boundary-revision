#include<bits/stdc++.h>
#include<unistd.h>
#define N 5000000
#define minIntronLength 100
#define MAXERROR 50

typedef long long ll;
using namespace std;

//**************设置各个参数*********************************************
int maxIntronFix;//intron小于这个值时合并两个模块
char* inputFilePath;
char* outputFilePath;//输入文件路径，输出文件路径
vector<char* > evdFilePath;//证据文件路径（可能有好几个，所以用vector）
//***********************************************************************



//***************定义PSL文件及BED文件的结构体*************************************
typedef struct PSLDefine
{
    ll equal_mRNA_to_DNA_number;			//0		mRNA与DNA完全相同碱基对数
    ll unequal_mRNA_to_DNA_number;		//1		mRNA与DNA不同碱基对数
    ll insert_mRNA_fragment;		//4		mRNA插入片段数
    ll insert_mRNA_bp_pair;		//5		(4)中片段碱基对总和
    ll deficiency_mRNA_fragment;	//6		mRNA相对于DNA缺失片段数
    ll deficiency_mRNA_bp_pair;	//7		(6)中片段碱基对总和
    string plus_minus_flag;		//8		正负链标志
    string mRNA_name;				//9		mRNA编号
    ll mRNA_bp_pair_number;		//10	mRNA碱基对总数
    ll mRNA_comparison_begin_index;	//11	相对于mRNA的开始，开始比对的位置(0-base)
    ll mRNA_comparison_end_index;		//12	相对于mRNA的开始，结束比对的位置(1-base)
    string DNA_name;					//13	DNA编号
    ll DNA_bp_pair_number;			//14	DNA碱基对总数
    ll DNA_comparison_begin_index;	//15	相对于DNA的开始，开始比对的位置(0-base)
    ll DNA_comparison_end_index;		//16	相对于DNA的开始，结束比对的位置(1-base)
    ll comparsion_block_number;				//17	模块总数
    vector<ll> block_bp_pair_number;		//18	模块的碱基对数(以逗号分隔)
    vector<ll> mRNA_block_begin_index;	//19	相对于mRNA的开始，模块开始的位置(0-base)
    vector<ll> DNA_block_begin_index;	//20	相对于DNA的开始，模块开始的位置(0-base)

} PSLDefine;

struct BEDDefine
{
    string DNA_name;					//唯一的DNA名字
    ll DNA_comparison_begin_index;
    ll DNA_comparison_end_index;
    string mRNA_name;					//mRNA名字由index|来源mRNA名字1|来源mRNA名字2|......
    ll same_mRNA_number;
    char plus_minus_flag;
    ll important_area_begin_index;
    ll important_area_end_index;
    char isempty;
    ll comparison_block_number;
    vector<ll> block_bp_pair_number;
    vector<ll> DNA_block_begin_index;
    bool operator<(const BEDDefine& a)const
    {
        if(DNA_name==a.DNA_name)
        {
            if(DNA_comparison_begin_index==a.DNA_comparison_begin_index)
            {
                if(DNA_comparison_end_index==a.DNA_comparison_end_index)
                {
                    return comparison_block_number<a.comparison_block_number;
                }
                return DNA_comparison_end_index<a.DNA_comparison_end_index;
            }
            return DNA_comparison_begin_index<a.DNA_comparison_begin_index;
        }
        return DNA_name<a.DNA_name;
    }
} ;
//******************************************************************************************



//***********************软件文件的结构体****************************************************
typedef struct Mystruct
{
    string name;//Acceptor or Dornor
    ll Dig1;
    ll Dig2;
    double score;//打分
    string strand;//符号

    bool operator < (const Mystruct a) const
    {
        return Dig1<a.Dig1;
    }
} DorAcc;
//*******************************************************************************************
//***********************************************
/*读四个类型文件*/  //属于软件的方法
set<Mystruct> AccZheng,AccFu,DorZheng,DorFu;



set<BEDDefine> result_fixIntron;//最后输出的bed文件

map<int,vector<string> > getmRNA_name;//第一步中相同的的mRNA的name

//***************一些函数***********************************
/*此函数将字符串以"_"分割放到字符串数组中*/
vector<string> name_Split(string m,char c)
{
    vector<string> kk;
    string a="";
    ll i;
    for(i=0;;)
    {
        while(1)
        {
            if(m[i]==c || i==m.size())
                break;
            a+=m[i];
            i++;

        }
        i++;
        kk.push_back(a);
        a="";
        if(i==m.size()+1)
            break;
    }
    return kk;
}

/*此函数是将字符串数组通过'_'连成一个字符串*/
string name_Combine(vector<string> kk)
{
    string m="";
    ll i;
    for(i=0; i<kk.size()-1; i++)
    {
        m+=kk[i];
        m+='_';
    }
    m+=kk[i];
    return m;
}
/*数字转成字符串*/
string Dig2Alp(ll a)
{
    /*stringstream m;//这种用法效率较低
    m<<a;
    return m.str();*/
    /*string str="";
    while(a){
        str+=a%10+'0';
        a/=10;
    }
    str=strrev(const_cast<char *>(str.c_str()));
    return str;*/
    string str="";
    string str_dig="";
    if(a==0)
    {
        str+="0";
    }
    while(a)
    {
        str_dig="";
        str_dig+=a%10+'0';
        str=str_dig+str;
        a/=10;
    }
    return str;


}
/*字符串转为数字*/
int Alp2Dig(string a)
{
    ll k;
    k=atoi(a.c_str());
    return k;
}

/*vector<string>转变为vector<int>*/
vector<ll> Vecs2Veci(vector<string> vecs)
{
    vector<ll> veci;
    vector<string>::iterator it;
    for(it=vecs.begin(); it!=vecs.end(); ++it)
        veci.push_back(Alp2Dig(*it));
    return veci;
}

//**********************************************************

//**************证据文件的结构体********************************
struct thr_Evidence
{
    ll s_low;//intron边界
    ll s_high;//intron边界
    string strance;//DNA正负链

    bool operator < (const thr_Evidence a) const
    {
        if(s_low==a.s_low)
            return s_high <a.s_high;
        return s_low<a.s_low;
    }
};

//**************************************************************


long long zyftest=0;

//**********************全局定义的参数***************************
ll all_index;
//定义两个数组保存证据
set<thr_Evidence> Evd_pos,Evd_nag;//正链，负链
//定义一个要修正的
//set<thr_Evidence> my_bed;
//定义一个map
//map<thr_Evidence,thr_Evidence> myMap;
struct Node
{
    ll start, b_end ;
    ll pos ;
    bool operator < (const Node a) const
    {
        if(start == a.start)
            return b_end > a.b_end ;
        return start < a.start ;
    }
};
BEDDefine bed_data_line;
vector<BEDDefine> bed_data;
int mmm=0;
PSLDefine in[N];
Node change[N] ;
map<string,ll> aa ;
ll bb[N] ;
ll ss[N] ;
ll cnt1 ;
ll cnt ;
ll te1[10000] ;
ll te2[10000] ;
bool vis[5000000];
map<string, vector<ll> > mode_one,mode_two;
//************************************************************************




//*******************读取证据文件*********************************
//****************************************************************
/*读证据的bed文件*/
void Evidence_bed(char* filepath)
{
    ifstream fin_ebed;
    string str;
    vector<string> vec_line;
    vec_line.clear();
    vector<ll> vec_block,vec_begin;//边界值，模块大小，模块开始位置开
    vec_block.clear();
    vec_begin.clear();
    ll st;
    thr_Evidence struct1;


    fin_ebed.open(filepath);
    //cout<<"begin of Evidence_bed()"<<endl;
    while(!fin_ebed.eof())
    {
        getline(fin_ebed,str);
        if(str=="")
            break;
        vec_line.clear();
        vec_line=name_Split(str,'\t');//得到一行的数组
        if(Alp2Dig(vec_line[9])==1)
            continue;
        else
        {
            st=Alp2Dig(vec_line[1]);

            vector<string> a;
            a.clear();
            a=name_Split(vec_line[10],',');
            a.pop_back();
            vec_block=Vecs2Veci(a);
            a.clear();
            a=name_Split(vec_line[11],',');
            a.pop_back();
            vec_begin=Vecs2Veci(a);

            //存储intron的边界
            for(int i=0; i<Alp2Dig(vec_line[9])-1; i++)
            {
                struct1.s_low=(vec_begin[i]+vec_block[i]+st);//边界
                struct1.s_high=(vec_begin[i+1]+st);//边界
                struct1.strance=vec_line[5];//正负号
                if(struct1.strance=="+")
                    Evd_pos.insert(struct1);
                else
                    Evd_nag.insert(struct1);
            }
        }
    }
    fin_ebed.close();
}

/*读证据的psl文件*/
void Evidence_psl(char* filepath)
{
    ifstream fin_epsl;
    string str;
    vector<string> vec_line;
    vec_line.clear();
    vector<ll> vec_block,vec_begin;//边界值，模块大小，模块开始位置开
    vec_block.clear();
    vec_begin.clear();
    ll st;
    thr_Evidence struct1;
    vector<ll> te1,te2;//记录分割后的字符
    te1.clear();
    te2.clear();

    fin_epsl.open(filepath);

    //fin_epsl.open("..\\2017.2.15\\mytest_evidence_psl.psl");
    //cout<<"begin of Evidence_psl()"<<endl;
    while(!fin_epsl.eof())
    {
        //cout<<mycount++<<endl;
        getline(fin_epsl,str);
        if(str=="")
            break;
        //cout<<zyftest++<<endl;
        //mycount+=1;
        vec_line.clear();
        vec_line=name_Split(str,'\t');
        if(Alp2Dig(vec_line[17])==1)
            continue;
        else
        {
            vector<string> a;
            a.clear();
            int ccc=0;
            a=name_Split(vec_line[18],',');
            a.pop_back();
            te1.clear();
            te1=Vecs2Veci(a);
            a.clear();
            a=name_Split(vec_line[20],',');
            a.pop_back();
            te2.clear();
            te2=Vecs2Veci(a);
            //先把能合并的进行合并
            vec_begin.clear();
            vec_block.clear();
            vec_block.push_back(te1[0]);
            vec_begin.push_back(te2[0]);
            for(int i=1; i<Alp2Dig(vec_line[17]); i++)
            {
                if(te2[i]<=vec_begin[ccc]+vec_block[ccc]+maxIntronFix)//最小值设定为9
                    vec_block[ccc]=te2[i]+te1[i]-vec_begin[ccc];
                else
                {
                    vec_block.push_back(te1[i]);
                    vec_begin.push_back(te2[i]);
                    ccc++;
                }
            }


            //模块数为ccc+1
            if(ccc+1<2)
                continue;
            else
            {
                for(int j=0; j<ccc; j++)
                {
                    struct1.s_low=vec_begin[j]+vec_block[j];//边界
                    struct1.s_high=vec_begin[j+1];//边界
                    struct1.strance=vec_line[8];//正负号

                    //cout<<struct1.s_low<<" "<<struct1.s_high<<" "<<struct1.strance<<endl;

                    if(struct1.strance=="+")
                        Evd_pos.insert(struct1);
                    else
                        Evd_nag.insert(struct1);
                }
            }
            vec_begin.clear();
            vec_block.clear();
            vec_line.clear();

        }
        //cout<<mycount++<<endl;
    }
    fin_epsl.close();
}

/*读证据的txt文件*/
void Evidence_txt(char* filepath)
{
    ifstream fin_etxt;
    string str;
    vector<string> vec_line;
    vec_line.clear();
    thr_Evidence struct1;


    fin_etxt.open(filepath);
    //cout<<"begin of Evidence_txt()"<<endl;
    while(!fin_etxt.eof())
    {
        getline(fin_etxt,str);
        if(str=="")
            break;
        vec_line.clear();
        vec_line=name_Split(str,'\t');
        if(name_Split(vec_line[4],'-').size()<4)
        {
            continue;
        }
        else
        {
            vector<string> vec_boundry;
            vec_boundry.clear();

            vec_boundry=name_Split(vec_line[4],'-');
            int kkk=1;
            for(int i=0; i<vec_boundry.size()/2-1; i++)
            {
                struct1.s_low=Alp2Dig(vec_boundry[kkk]);
                struct1.s_high=Alp2Dig(vec_boundry[kkk+1]);
                struct1.strance=vec_line[5];
                if(struct1.strance=="+")
                    Evd_pos.insert(struct1);
                else
                    Evd_nag.insert(struct1);
                kkk+=2;
            }
        }
    }
    fin_etxt.close();
}

//****************************************************************
