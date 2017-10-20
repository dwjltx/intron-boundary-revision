#include<bits/stdc++.h>
#include<unistd.h>
#define N 5000000
#define minIntronLength 100
#define MAXERROR 50

typedef long long ll;
using namespace std;

//**************���ø�������*********************************************
int maxIntronFix;//intronС�����ֵʱ�ϲ�����ģ��
char* inputFilePath;
char* outputFilePath;//�����ļ�·��������ļ�·��
vector<char* > evdFilePath;//֤���ļ�·���������кü�����������vector��
//***********************************************************************



//***************����PSL�ļ���BED�ļ��Ľṹ��*************************************
typedef struct PSLDefine
{
    ll equal_mRNA_to_DNA_number;			//0		mRNA��DNA��ȫ��ͬ�������
    ll unequal_mRNA_to_DNA_number;		//1		mRNA��DNA��ͬ�������
    ll insert_mRNA_fragment;		//4		mRNA����Ƭ����
    ll insert_mRNA_bp_pair;		//5		(4)��Ƭ�μ�����ܺ�
    ll deficiency_mRNA_fragment;	//6		mRNA�����DNAȱʧƬ����
    ll deficiency_mRNA_bp_pair;	//7		(6)��Ƭ�μ�����ܺ�
    string plus_minus_flag;		//8		��������־
    string mRNA_name;				//9		mRNA���
    ll mRNA_bp_pair_number;		//10	mRNA���������
    ll mRNA_comparison_begin_index;	//11	�����mRNA�Ŀ�ʼ����ʼ�ȶԵ�λ��(0-base)
    ll mRNA_comparison_end_index;		//12	�����mRNA�Ŀ�ʼ�������ȶԵ�λ��(1-base)
    string DNA_name;					//13	DNA���
    ll DNA_bp_pair_number;			//14	DNA���������
    ll DNA_comparison_begin_index;	//15	�����DNA�Ŀ�ʼ����ʼ�ȶԵ�λ��(0-base)
    ll DNA_comparison_end_index;		//16	�����DNA�Ŀ�ʼ�������ȶԵ�λ��(1-base)
    ll comparsion_block_number;				//17	ģ������
    vector<ll> block_bp_pair_number;		//18	ģ��ļ������(�Զ��ŷָ�)
    vector<ll> mRNA_block_begin_index;	//19	�����mRNA�Ŀ�ʼ��ģ�鿪ʼ��λ��(0-base)
    vector<ll> DNA_block_begin_index;	//20	�����DNA�Ŀ�ʼ��ģ�鿪ʼ��λ��(0-base)

} PSLDefine;

struct BEDDefine
{
    string DNA_name;					//Ψһ��DNA����
    ll DNA_comparison_begin_index;
    ll DNA_comparison_end_index;
    string mRNA_name;					//mRNA������index|��ԴmRNA����1|��ԴmRNA����2|......
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



//***********************����ļ��Ľṹ��****************************************************
typedef struct Mystruct
{
    string name;//Acceptor or Dornor
    ll Dig1;
    ll Dig2;
    double score;//���
    string strand;//����

    bool operator < (const Mystruct a) const
    {
        return Dig1<a.Dig1;
    }
} DorAcc;
//*******************************************************************************************
//***********************************************
/*���ĸ������ļ�*/  //��������ķ���
set<Mystruct> AccZheng,AccFu,DorZheng,DorFu;



set<BEDDefine> result_fixIntron;//��������bed�ļ�

map<int,vector<string> > getmRNA_name;//��һ������ͬ�ĵ�mRNA��name

//***************һЩ����***********************************
/*�˺������ַ�����"_"�ָ�ŵ��ַ���������*/
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

/*�˺����ǽ��ַ�������ͨ��'_'����һ���ַ���*/
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
/*����ת���ַ���*/
string Dig2Alp(ll a)
{
    /*stringstream m;//�����÷�Ч�ʽϵ�
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
/*�ַ���תΪ����*/
int Alp2Dig(string a)
{
    ll k;
    k=atoi(a.c_str());
    return k;
}

/*vector<string>ת��Ϊvector<int>*/
vector<ll> Vecs2Veci(vector<string> vecs)
{
    vector<ll> veci;
    vector<string>::iterator it;
    for(it=vecs.begin(); it!=vecs.end(); ++it)
        veci.push_back(Alp2Dig(*it));
    return veci;
}

//**********************************************************

//**************֤���ļ��Ľṹ��********************************
struct thr_Evidence
{
    ll s_low;//intron�߽�
    ll s_high;//intron�߽�
    string strance;//DNA������

    bool operator < (const thr_Evidence a) const
    {
        if(s_low==a.s_low)
            return s_high <a.s_high;
        return s_low<a.s_low;
    }
};

//**************************************************************


long long zyftest=0;

//**********************ȫ�ֶ���Ĳ���***************************
ll all_index;
//�����������鱣��֤��
set<thr_Evidence> Evd_pos,Evd_nag;//����������
//����һ��Ҫ������
//set<thr_Evidence> my_bed;
//����һ��map
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




//*******************��ȡ֤���ļ�*********************************
//****************************************************************
/*��֤�ݵ�bed�ļ�*/
void Evidence_bed(char* filepath)
{
    ifstream fin_ebed;
    string str;
    vector<string> vec_line;
    vec_line.clear();
    vector<ll> vec_block,vec_begin;//�߽�ֵ��ģ���С��ģ�鿪ʼλ�ÿ�
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
        vec_line=name_Split(str,'\t');//�õ�һ�е�����
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

            //�洢intron�ı߽�
            for(int i=0; i<Alp2Dig(vec_line[9])-1; i++)
            {
                struct1.s_low=(vec_begin[i]+vec_block[i]+st);//�߽�
                struct1.s_high=(vec_begin[i+1]+st);//�߽�
                struct1.strance=vec_line[5];//������
                if(struct1.strance=="+")
                    Evd_pos.insert(struct1);
                else
                    Evd_nag.insert(struct1);
            }
        }
    }
    fin_ebed.close();
}

/*��֤�ݵ�psl�ļ�*/
void Evidence_psl(char* filepath)
{
    ifstream fin_epsl;
    string str;
    vector<string> vec_line;
    vec_line.clear();
    vector<ll> vec_block,vec_begin;//�߽�ֵ��ģ���С��ģ�鿪ʼλ�ÿ�
    vec_block.clear();
    vec_begin.clear();
    ll st;
    thr_Evidence struct1;
    vector<ll> te1,te2;//��¼�ָ����ַ�
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
            //�Ȱ��ܺϲ��Ľ��кϲ�
            vec_begin.clear();
            vec_block.clear();
            vec_block.push_back(te1[0]);
            vec_begin.push_back(te2[0]);
            for(int i=1; i<Alp2Dig(vec_line[17]); i++)
            {
                if(te2[i]<=vec_begin[ccc]+vec_block[ccc]+maxIntronFix)//��Сֵ�趨Ϊ9
                    vec_block[ccc]=te2[i]+te1[i]-vec_begin[ccc];
                else
                {
                    vec_block.push_back(te1[i]);
                    vec_begin.push_back(te2[i]);
                    ccc++;
                }
            }


            //ģ����Ϊccc+1
            if(ccc+1<2)
                continue;
            else
            {
                for(int j=0; j<ccc; j++)
                {
                    struct1.s_low=vec_begin[j]+vec_block[j];//�߽�
                    struct1.s_high=vec_begin[j+1];//�߽�
                    struct1.strance=vec_line[8];//������

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

/*��֤�ݵ�txt�ļ�*/
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
