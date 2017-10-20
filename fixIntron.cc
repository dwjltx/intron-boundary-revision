#include"fixIntron.h"

//读取输入路径的文件
void Read()
{


    PSLDefine temp ;
    ll rubbish0, rubbish1 ;
    mode_one.clear();
    mode_two.clear();
    cnt = 0 ;
    string hh="\t";
    char cc;
    while(scanf("%lld\t", &in[cnt].equal_mRNA_to_DNA_number) != EOF)
    {
        scanf("%lld\t", &in[cnt].unequal_mRNA_to_DNA_number) ;
        scanf("%lld\t%lld\t", &rubbish0, &rubbish1) ;
        scanf("%lld\t%lld\t", &in[cnt].insert_mRNA_fragment, &in[cnt].insert_mRNA_bp_pair);
        scanf("%lld\t%lld\t", &in[cnt].deficiency_mRNA_fragment, &in[cnt].deficiency_mRNA_bp_pair) ;
        cin >> in[cnt].plus_minus_flag;
        cc=getchar();
        cin>> in[cnt].mRNA_name ;
        cc=getchar();
        scanf("%lld\t%lld\t%lld\t", &in[cnt].mRNA_bp_pair_number, &in[cnt].mRNA_comparison_begin_index, &in[cnt].mRNA_comparison_end_index) ;
        cin >> in[cnt].DNA_name ;
        cc=getchar();
        scanf("%lld\t%lld\t%lld\t", &in[cnt].DNA_bp_pair_number, &in[cnt].DNA_comparison_begin_index, &in[cnt].DNA_comparison_end_index) ;
        scanf("%lld\t", &in[cnt].comparsion_block_number) ;
        for(ll i = 0 ; i < in[cnt].comparsion_block_number ; i ++)
        {
            ll tt ;
            scanf("%lld,",&tt) ;
            te1[i] = tt ;
        }
        cc=getchar();

        for(ll i = 0 ; i < in[cnt].comparsion_block_number ; i ++)
        {
            ll tt ;
            scanf("%lld,",&tt) ;
            in[cnt].mRNA_block_begin_index.push_back(tt) ;
        }
        cc=getchar();
        for(ll i = 0 ; i < in[cnt].comparsion_block_number ; i ++)
        {
            ll tt ;
            scanf("%lld,",&tt) ;
            te2[i] = tt ;
        }
        cc=getchar();
        ll ccc = 0 ;
        in[cnt].DNA_block_begin_index.push_back(te2[0]) ;
        in[cnt].block_bp_pair_number.push_back(te1[0]) ;
        //合并所需
        for(ll i = 1 ; i < in[cnt].comparsion_block_number ; i ++)
        {
            if(te2[i] <= in[cnt].DNA_block_begin_index[ccc] + in[cnt].block_bp_pair_number[ccc] + maxIntronFix)
            {
                in[cnt].block_bp_pair_number[ccc] = te2[i] + te1[i] - in[cnt].DNA_block_begin_index[ccc] ;
            }
            else
            {
                in[cnt].DNA_block_begin_index.push_back(te2[i]) ;
                in[cnt].block_bp_pair_number.push_back(te1[i]) ;
                ccc ++ ;
            }
        }

        in[cnt].comparsion_block_number = ccc + 1;

        string my_string=in[cnt].plus_minus_flag+in[cnt].DNA_name;
        if(in[cnt].comparsion_block_number==1)
            mode_one[my_string].push_back(cnt);
        if(in[cnt].comparsion_block_number>1)
            mode_two[my_string].push_back(cnt);
        cnt ++ ;
    }
}


//一个模块的
void solve_one(vector<ll> one)
{
    int l=one.size();
    memset(vis,false,sizeof(vis));
    if(l!=0)
    {
        for(int i=0; i<l; i++)
        {
            int t=one[i];
            change[i].start=in[t].DNA_block_begin_index[0];
            change[i].b_end=in[t].DNA_block_begin_index[0]+in[t].block_bp_pair_number[0];
            change[i].pos=t;
        }
        sort(change,change+l);
        ll num = 1 ;
        ll left = change[0].start ;
        ll right = change[0].b_end ;
        for(ll i = 1 ; i < l ; i ++)
        {
            if(change[i].start > right)
            {
                getmRNA_name[all_index].push_back(in[one[i-1] ].mRNA_name);

                bed_data_line.DNA_name=in[change[i-1].pos].DNA_name;
                bed_data_line.DNA_comparison_begin_index=left;
                bed_data_line.DNA_comparison_end_index=right;
                bed_data_line.mRNA_name="ioseq"+Dig2Alp(all_index);
                //bed_data_line.mRNA_name+="|";
                //bed_data_line.mRNA_name+=in[one[i-1] ].mRNA_name;
                if(getmRNA_name[all_index].size()>20)
                {
                    for(int im=0; im<20; im++)
                    {
                        bed_data_line.mRNA_name+="|";
                        bed_data_line.mRNA_name+=getmRNA_name[all_index][im];
                    }
                }
                else
                {
                    for(int im=0; im<getmRNA_name[all_index].size(); im++)
                    {
                        bed_data_line.mRNA_name+="|";
                        bed_data_line.mRNA_name+=getmRNA_name[all_index][im];
                    }
                }


                bed_data_line.same_mRNA_number=num;
                if(in[one[0]].plus_minus_flag=="+")
                    bed_data_line.plus_minus_flag='+';
                else
                    bed_data_line.plus_minus_flag='-';
                bed_data_line.important_area_begin_index=left;
                bed_data_line.important_area_end_index=right;
                bed_data_line.isempty='.';
                bed_data_line.comparison_block_number=1;
                bed_data_line.block_bp_pair_number.push_back(right-left);
                bed_data_line.DNA_block_begin_index.push_back(0);
                bed_data.push_back(bed_data_line);//保存起来
                bed_data_line.block_bp_pair_number.clear();
                bed_data_line.DNA_block_begin_index.clear();

                all_index++;
                num = 1 ;
                left = change[i].start ;
                right = change[i].b_end ;
            }
            else
            {
                left = min(left,change[i].start) ;
                right = max(right,change[i].b_end) ;
                num ++ ;
                getmRNA_name[all_index].push_back(in[one[i-1] ].mRNA_name);
            }
        }

        getmRNA_name[all_index].push_back(in[one[l-1] ].mRNA_name);
        bed_data_line.DNA_name=in[change[l-1].pos].DNA_name;
        bed_data_line.DNA_comparison_begin_index=left;
        bed_data_line.DNA_comparison_end_index=right;
        bed_data_line.mRNA_name="ioseq"+Dig2Alp(all_index);
        //bed_data_line.mRNA_name+="|";
        //bed_data_line.mRNA_name+=in[one[0] ].mRNA_name;
        if(getmRNA_name[all_index].size()>20)
        {
            for(int im=0; im<20; im++)
            {
                bed_data_line.mRNA_name+="|";
                bed_data_line.mRNA_name+=getmRNA_name[all_index][im];
            }
        }
        else
        {
            for(int im=0; im<getmRNA_name[all_index].size(); im++)
            {
                bed_data_line.mRNA_name+="|";
                bed_data_line.mRNA_name+=getmRNA_name[all_index][im];
            }
        }


        bed_data_line.same_mRNA_number=num;
        if(in[one[0]].plus_minus_flag=="+")
            bed_data_line.plus_minus_flag='+';
        else
            bed_data_line.plus_minus_flag='-';
        bed_data_line.important_area_begin_index=left;
        bed_data_line.important_area_end_index=right;
        bed_data_line.isempty='.';
        bed_data_line.comparison_block_number=1;
        bed_data_line.block_bp_pair_number.push_back(right-left);
        bed_data_line.DNA_block_begin_index.push_back(0);
        bed_data.push_back(bed_data_line);//保存起来
        bed_data_line.block_bp_pair_number.clear();
        bed_data_line.DNA_block_begin_index.clear();
        all_index++;
    }
}

//两个模块的
void solve_two(vector<ll> two)
{
    getmRNA_name.clear();
    memset(bb,0,sizeof bb) ;
    memset(ss,0,sizeof ss) ;
    aa.clear() ;
    cnt1 = 1 ;
    ll l = two.size() ;
    for(ll i = 0 ; i < l ; i ++)
    {
        ll t = two[i] ;
        ll left = in[t].DNA_block_begin_index[0] ;
        ll num = in[t].DNA_block_begin_index.size() ;
        ll right = in[t].DNA_block_begin_index[num-1] + in[t].block_bp_pair_number[num-1] ;
        string ttt = "" ;
        ll num_temp = (left+in[t].block_bp_pair_number[0]) ;
        while(num_temp)
        {
            ttt += num_temp%10+'0' ;
            num_temp/=10 ;
        }
        ttt += ',' ;
        for(ll j = 1 ; j < in[t].comparsion_block_number - 1 ; j ++)
        {
            num_temp = in[t].DNA_block_begin_index[j] ;
            while(num_temp)
            {
                ttt += num_temp%10+'0' ;
                num_temp/=10 ;
            }
            ttt += ',' ;
            num_temp = in[t].block_bp_pair_number[j] ;
            while(num_temp)
            {
                ttt += num_temp%10+'0' ;
                num_temp/=10 ;
            }
            ttt += ',' ;
        }
        num_temp = in[t].DNA_block_begin_index[num-1] ;
        while(num_temp)
        {
            ttt += num_temp%10+'0' ;
            num_temp/=10 ;
        }
        ttt += ',' ;
        if(aa[ttt] == 0)
        {
            getmRNA_name[cnt1].push_back(in[t].mRNA_name);//得到mRNA_name
            aa[ttt] = cnt1 ;
            bb[cnt1] ++ ;
            ss[cnt1] = t ;
            change[cnt1].start = left ;
            change[cnt1].b_end = right ;
            cnt1 ++ ;
        }
        else
        {

            ll temp = aa[ttt] ;
            getmRNA_name[temp].push_back(in[t].mRNA_name);//得到mRNA_name
            bb[temp] ++ ;
            change[temp].start = min(change[temp].start,left) ;
            change[temp].b_end = max(change[temp].b_end,right) ;
        }
    }
    for(ll i = 1 ; i < cnt1 ; i ++)
    {
        ll tt = in[ss[i]].comparsion_block_number ;
        bed_data_line.DNA_name=in[ss[i]].DNA_name;
        bed_data_line.DNA_comparison_begin_index=change[i].start;
        bed_data_line.DNA_comparison_end_index=change[i].b_end;
        bed_data_line.mRNA_name="ioseq"+Dig2Alp(all_index);
        for(int im=0; im<getmRNA_name[i].size(); im++)
        {
            bed_data_line.mRNA_name+="|";
            bed_data_line.mRNA_name+=getmRNA_name[i][im];
        }

        bed_data_line.same_mRNA_number=bb[i];
        if(in[two[0]].plus_minus_flag=="+")
            bed_data_line.plus_minus_flag='+';
        else
            bed_data_line.plus_minus_flag='-';
        bed_data_line.important_area_begin_index=change[i].start;
        bed_data_line.important_area_end_index=change[i].b_end;
        bed_data_line.isempty='.';
        bed_data_line.comparison_block_number=tt;
        all_index++;
        bed_data_line.block_bp_pair_number.push_back(in[ss[i]].block_bp_pair_number[0] + in[ss[i]].DNA_block_begin_index[0] - change[i].start);
        for(ll j = 1 ; j < tt - 1 ; j ++)
        {
            bed_data_line.block_bp_pair_number.push_back(in[ss[i]].block_bp_pair_number[j]);
        }
        bed_data_line.block_bp_pair_number.push_back(change[i].b_end-in[ss[i]].DNA_block_begin_index[tt-1]);
        bed_data_line.DNA_block_begin_index.push_back(change[i].start - change[i].start);
        for(ll j = 1 ; j < tt ; j ++)
        {
            bed_data_line.DNA_block_begin_index.push_back(in[ss[i]].DNA_block_begin_index[j]-change[i].start);
        }
        bed_data.push_back(bed_data_line);
        bed_data_line.block_bp_pair_number.clear();
        bed_data_line.DNA_block_begin_index.clear();
    }
}

void Solve()
{
    map<string,vector<ll> >::iterator it;
    for(it=mode_one.begin(); it!=mode_one.end(); ++it)
        solve_one(it->second);
    getmRNA_name.clear();
    for(it=mode_two.begin(); it!=mode_two.end(); ++it)
        solve_two(it->second);
    getmRNA_name.clear();
}

//bed文件写成文件
void Print_bed(BEDDefine bed_line,ofstream& out)
{
    out<<bed_line.DNA_name<<"\t";
    out<<bed_line.DNA_comparison_begin_index<<"\t"<<bed_line.DNA_comparison_end_index<<"\t";
    out<<bed_line.mRNA_name<<"\t";
    out<<bed_line.same_mRNA_number<<"\t";
    out<<bed_line.plus_minus_flag<<"\t";
    out<<bed_line.important_area_begin_index<<"\t"<<bed_line.important_area_end_index<<"\t";
    out<<".\t"<<bed_line.comparison_block_number<<"\t";
    for(ll j=0; j<bed_line.comparison_block_number; j++)
    {
        out<<bed_line.block_bp_pair_number[j]<<",";
    }
    out<<"\t";
    for(ll j=0; j<bed_line.comparison_block_number; j++)
    {
        out<<bed_line.DNA_block_begin_index[j]<<",";
    }
    out<<"\n";
}

//*****************读取所有所有证据文件********************************
void ReadAllEvd(vector<char* > evdFile)
{
    string filepath;
    //cout<<"beging function of ReadAllEvd()"<<endl;
    for(int i=0; i<evdFile.size(); i++)
    {
        filepath=evdFile[i];
        if(filepath[filepath.size()-1]=='d') //bed文件
        {
            Evidence_bed(evdFile[i]);
            //cout<<"deal bed evidence"<<endl;
        }
        else if(filepath[filepath.size()-1]=='l') //psl文件
        {
            Evidence_psl(evdFile[i]);
            //cout<<"deal psl evidence"<<endl;
        }
        else if(filepath[filepath.size()-1]=='t') //txt文件
        {
            Evidence_txt(evdFile[i]);
            //cout<<"deal txt evidence"<<endl;
        }
        else
        {
            cout<<"证据文件错误"<<endl;
            //文件类型错误
        }
    }
}
//**********************************************************************

//有证据文件的
void SecondStep()
{
    ofstream out(outputFilePath);
    BEDDefine bed_line,my_line;
    string str;
    int n;
    result_fixIntron.clear();
    //ReadEvidence();
    //cout<<"bbbbbbbbbbbbbb"<<endl;
    //cout<<evdFilePath.size()<<endl;
    ReadAllEvd(evdFilePath);
    //cout<<"aaaaaaaaaaaaaa"<<endl;
    //cout<<"Evd_pos.size::"<<Evd_pos.size()<<endl;
    vector<int> fixed_intron;//记录修正后的值
    int intron_begin,intron_end;//记录内含子的边界
    int intron_begin_after,intron_end_after;//记录前后修正后的值
    //cout<<"aaaaaaaaaaaaa"<<endl;
    //ReadEvidence();
    //cout<<"bed_data.size::"<<bed_data.size()<<endl;

    for(int i1=0; i1<bed_data.size(); i1++)
    {

        bed_line=bed_data[i1];
        fixed_intron.clear();
        //zyf++;
        //my_line=bed_line;
        // cout<<++zyf<<endl;
        // }
        //cout<<"张宜飞"<<bed_line.comparison_block_number<<endl;
        //*************输出结构体测试******************
        //PrintBED(bed_line);
        //*********************************************
        //PrintBED(bed_line);
        if(bed_line.comparison_block_number<2)
        {
            Print_bed(bed_line,out);
        }//当只有一个模块时直接输出，不用处理
        //PrintBED(bed_line);
        else
        {
            int i;
            string mystrance;//正负链
            int sum_low_high=888888;
            thr_Evidence First_Struct,Median_Struct;
            for(i=0; i<bed_line.comparison_block_number-1; i++) //循环的次数即内含子个数
            {
                intron_begin=bed_line.DNA_block_begin_index[i]+bed_line.block_bp_pair_number[i]+bed_line.DNA_comparison_begin_index;
                intron_end=bed_line.DNA_block_begin_index[i+1]+bed_line.DNA_comparison_begin_index;
                //输出内含子边界测试
                //cout<<intron_begin<<'\t'<<intron_end<<endl;
                if(bed_line.plus_minus_flag=='+')
                {
                    Median_Struct.s_low=intron_begin;
                    Median_Struct.s_high=intron_end;
                    Median_Struct.strance="+";

                    set<thr_Evidence>::iterator zhit,zhit_jia,zhit_jian;
                    First_Struct.s_low=intron_begin;
                    First_Struct.s_high=intron_end;
                    First_Struct.strance="+";

                    zhit=Evd_pos.lower_bound(First_Struct);//找到最近边界的那个值的迭代器
                    zhit_jia=zhit;
                    zhit_jian=zhit;
                    while(1)   //向上加的
                    {
                        if(abs((*zhit_jia).s_low-intron_begin)>minIntronLength)
                            break;
                        if(abs((*zhit_jia).s_high-intron_end)<minIntronLength && abs((*zhit_jia).s_low-intron_begin)+abs((*zhit_jia).s_high-intron_end)<sum_low_high)
                        {
                            sum_low_high=abs((*zhit_jia).s_low-intron_begin)+abs((*zhit_jia).s_high-intron_end);
                            Median_Struct=(*zhit_jia);
                        }
                        zhit_jia++;
                    }
                    while(1)   //向下减的
                    {
                        if(abs((*zhit_jian).s_low-intron_begin)>minIntronLength)
                            break;
                        if(abs((*zhit_jian).s_high-intron_end)<minIntronLength && abs((*zhit_jian).s_low-intron_begin)+abs((*zhit_jian).s_high-intron_end)<sum_low_high)
                        {
                            sum_low_high=abs((*zhit_jian).s_low-intron_begin)+abs((*zhit_jian).s_high-intron_end);
                            Median_Struct=(*zhit_jian);
                        }
                        zhit_jian--;
                    }
                    //至此，附件的遍历结束,最小值寻在中间结构体中
                    fixed_intron.push_back(Median_Struct.s_low);
                    fixed_intron.push_back(Median_Struct.s_high);
                    sum_low_high=888888;
                }

                else if(bed_line.plus_minus_flag=='-')
                {
                    Median_Struct.s_low=intron_begin;
                    Median_Struct.s_high=intron_end;
                    Median_Struct.strance="-";

                    set<thr_Evidence>::iterator fuit,fuit_jia,fuit_jian;
                    First_Struct.s_low=intron_begin;
                    First_Struct.s_high=intron_end;
                    First_Struct.strance="-";

                    fuit=Evd_nag.lower_bound(First_Struct);//找到最近边界的那个值的迭代器
                    fuit_jia=fuit;
                    fuit_jian=fuit;
                    while(1)   //向上加的
                    {
                        if(abs((*fuit_jia).s_low-intron_begin)>minIntronLength)
                            break;
                        if(abs((*fuit_jia).s_high-intron_end)<minIntronLength && abs((*fuit_jia).s_low-intron_begin)+abs((*fuit_jia).s_high-intron_end)<sum_low_high)
                        {
                            sum_low_high=abs((*fuit_jia).s_low-intron_begin)+abs((*fuit_jia).s_high-intron_end);
                            Median_Struct=(*fuit_jia);
                        }
                        fuit_jia++;
                    }
                    while(1)   //向下减的
                    {
                        if(abs((*fuit_jian).s_low-intron_begin)>minIntronLength)
                            break;
                        if(abs((*fuit_jian).s_high-intron_end)<minIntronLength && abs((*fuit_jian).s_low-intron_begin)+abs((*fuit_jian).s_high-intron_end)<sum_low_high)
                        {
                            sum_low_high=abs((*fuit_jian).s_low-intron_begin)+abs((*fuit_jian).s_high-intron_end);
                            Median_Struct=(*fuit_jian);
                        }
                        fuit_jian--;
                    }
                    //至此，附件的遍历结束,最小值寻在中间结构体中
                    fixed_intron.push_back(Median_Struct.s_low);
                    fixed_intron.push_back(Median_Struct.s_high);
                    sum_low_high=888888;
                }

                else
                    cout<<"Error in plus_minus_flag"<<endl;

            }
            //修正完一行
            bed_line.block_bp_pair_number.clear();
            bed_line.DNA_block_begin_index.clear();
            int k=1;//k来索引new_intron
            bed_line.DNA_block_begin_index.push_back(0);
            bed_line.block_bp_pair_number.push_back(fixed_intron[0]-bed_line.DNA_comparison_begin_index);
            for(i=0; i<bed_line.comparison_block_number-2; i++)
            {
                bed_line.DNA_block_begin_index.push_back(fixed_intron[k]-bed_line.DNA_comparison_begin_index);
                bed_line.block_bp_pair_number.push_back(fixed_intron[k+1]-fixed_intron[k]);
                k+=2;
            }
            bed_line.DNA_block_begin_index.push_back(fixed_intron[k]-bed_line.DNA_comparison_begin_index);
            bed_line.block_bp_pair_number.push_back(bed_line.DNA_comparison_end_index-fixed_intron[k]);

            //把需要合并的合并

            vector<ll> tt1,tt2;//记录模块大小及开始位置
            tt1.clear();
            tt2.clear();
            tt1=bed_line.block_bp_pair_number;
            tt2=bed_line.DNA_block_begin_index;
            bed_line.block_bp_pair_number.clear();
            bed_line.DNA_block_begin_index.clear();
            int cc=0;//记录模块数
            bed_line.block_bp_pair_number.push_back(tt1[0]);
            bed_line.DNA_block_begin_index.push_back(tt2[0]);
            for(i=0; i<bed_line.comparison_block_number; i++)
            {
                if(tt2[i]<=bed_line.DNA_block_begin_index[cc]+bed_line.block_bp_pair_number[cc]+maxIntronFix)
                    bed_line.block_bp_pair_number[cc]=tt2[i]+tt1[i]-bed_line.DNA_block_begin_index[cc];
                else
                {
                    bed_line.block_bp_pair_number.push_back(tt1[i]);
                    bed_line.DNA_block_begin_index.push_back(tt2[i]);
                    cc++;
                }
            }
            bed_line.comparison_block_number=cc+1;
            //Print_bed(bed_line,out);//原先是直接输出，现在最后的文件需要排序
            result_fixIntron.insert(bed_line);
        }
    }
    //输出排序、去重之后的文件s数据
    set<BEDDefine>::iterator fixed_it;
    for(fixed_it=result_fixIntron.begin(); fixed_it!=result_fixIntron.end(); fixed_it++)
    {
        Print_bed(*fixed_it,out);
    }

}


//软件的
//读四个类型，用软件通过.fa文件得到的结果
void ReadAccDor()
{
    FILE *aczh=fopen("acceptor_zheng.txt","r");
    FILE *acfu=fopen("acceptor_fu.txt","r");
    FILE *dozh=fopen("dornor_zheng.txt","r");
    FILE *dofu=fopen("dornor_fu.txt","r");
    Mystruct str;
    char s1[10],s2[10],s3[10],s4[10],s5[10];//结构体内部
    while(fscanf(aczh,"%s\t%s\t%s\t%s\t%s",s1,s2,s3,s4,s5)!=EOF){
        str.name=(string)s1;
        str.Dig1=Alp2Dig((string)s2);
        str.Dig2=Alp2Dig((string)s3);
        str.score=atof(s4);
        str.strand=(string)s5;
        //PrintStru(str);
        AccZheng.insert(str);
    }
    while(fscanf(acfu,"%s\t%s\t%s\t%s\t%s",s1,s2,s3,s4,s5)!=EOF){
        str.name=(string)s1;
        str.Dig1=Alp2Dig((string)s2);
        str.Dig2=Alp2Dig((string)s3);
        str.score=atof(s4);
        str.strand=(string)s5;
        //PrintStru(str);
        AccFu.insert(str);
    }
    while(fscanf(dozh,"%s\t%s\t%s\t%s\t%s",s1,s2,s3,s4,s5)!=EOF){
        str.name=(string)s1;
        str.Dig1=Alp2Dig((string)s2);
        str.Dig2=Alp2Dig((string)s3);
        str.score=atof(s4);
        str.strand=(string)s5;
        //PrintStru(str);
        DorZheng.insert(str);
    }
    while(fscanf(dofu,"%s\t%s\t%s\t%s\t%s",s1,s2,s3,s4,s5)!=EOF){
        str.name=(string)s1;
        str.Dig1=Alp2Dig((string)s2);
        str.Dig2=Alp2Dig((string)s3);
        str.score=atof(s4);
        str.strand=(string)s5;
        //PrintStru(str);
        DorFu.insert(str);
    }
}

//没有证据文件的，软件得到的对比
void SecondStep_2()
{

    ReadAccDor();//读完之后得到四类证据，约30s
    ofstream out(outputFilePath);
    BEDDefine bed_line;
    string str;
    int n;

    vector<int> fixed_intron;//记录修正后的值
    int intron_begin,intron_end;//记录内含子的边界
    int intron_begin_after,intron_end_after;//记录前后修正后的值
    for(ll i1=0;i1<bed_data.size();i1++){
            bed_line=bed_data[i1];
            fixed_intron.clear();


            if(bed_line.comparison_block_number<2){
                Print_bed(bed_line,out);
            }//当只有一个模块时直接输出，不用处理
                //PrintBED(bed_line);
            else{
                int i;
                string mystrance;//正负链
                int sum_low_high=8888;
                Mystruct MedStru;//中间过程的结构体
                for(i=0;i<bed_line.comparison_block_number-1;i++){//循环的次数即内含子个数
                    intron_begin=bed_line.DNA_block_begin_index[i]+bed_line.block_bp_pair_number[i]+bed_line.DNA_comparison_begin_index;
                    intron_end=bed_line.DNA_block_begin_index[i+1]+bed_line.DNA_comparison_begin_index;
                    //输出内含子边界测试
                    //cout<<intron_begin<<'\t'<<intron_end<<endl;

                    if(bed_line.plus_minus_flag=='+')//当是正链时
                    {
                        //正链前面的是Donor
                        MedStru.name="Donor";
                        MedStru.Dig1=intron_begin;
                        MedStru.Dig2=intron_begin+1;
                        MedStru.score=0.0;
                        MedStru.strand="+";

                        set<Mystruct>::iterator it,it1;
                        it=DorZheng.lower_bound(MedStru);
                        it1=it;
                        it1--;
                        if(abs((*it1).Dig1-intron_begin)>MAXERROR && abs((*it).Dig1-intron_begin)>MAXERROR)
                            fixed_intron.push_back(intron_begin);
                        else{
                            if(abs((*it1).Dig1-intron_begin)<=abs((*it).Dig1-intron_begin))
                                fixed_intron.push_back((*it1).Dig1);
                            else
                                fixed_intron.push_back((*it).Dig1);
                        }
                        //正链后面的是Acceptor
                        MedStru.name="Acceptor";
                        MedStru.Dig1=intron_end-1;
                        MedStru.Dig2=intron_end;
                        MedStru.score=0.0;
                        MedStru.strand="+";
                        it=AccZheng.lower_bound(MedStru);
                        it1=it;
                        it1--;
                        if(abs((*it1).Dig1-intron_end)>MAXERROR && abs((*it).Dig1-intron_end)>MAXERROR)
                            fixed_intron.push_back(intron_end);
                        else{
                            if(abs((*it1).Dig1-intron_end)<=abs((*it).Dig1-intron_end))
                                fixed_intron.push_back((*it1).Dig1);
                            else
                                fixed_intron.push_back((*it).Dig1);
                        }


                    }
                    else if(bed_line.plus_minus_flag=='-')//当是负链时
                    {
                       //负链前面的是Acceptor
                        MedStru.name="Acceptor";
                        MedStru.Dig1=intron_begin;
                        MedStru.Dig2=intron_begin+1;
                        MedStru.score=0.0;
                        MedStru.strand="-";

                        set<Mystruct>::iterator it,it1;
                        it=DorZheng.lower_bound(MedStru);
                        it1=it;
                        it1--;
                        if(abs((*it1).Dig1-intron_begin)>MAXERROR && abs((*it).Dig1-intron_begin)>MAXERROR)
                            fixed_intron.push_back(intron_begin);
                        else{
                            if(abs((*it1).Dig1-intron_begin)<=abs((*it).Dig1-intron_begin))
                                fixed_intron.push_back((*it1).Dig1);
                            else
                                fixed_intron.push_back((*it).Dig1);
                        }
                        //负链后面的是Donor
                        MedStru.name="Donor";
                        MedStru.Dig1=intron_end-1;
                        MedStru.Dig2=intron_end;
                        MedStru.score=0.0;
                        MedStru.strand="-";
                        it=AccZheng.lower_bound(MedStru);
                        it1=it;
                        it1--;
                        if(abs((*it1).Dig1-intron_end)>MAXERROR && abs((*it).Dig1-intron_end)>MAXERROR)
                            fixed_intron.push_back(intron_end);
                        else{
                            if(abs((*it1).Dig1-intron_end)<=abs((*it).Dig1-intron_end))
                                fixed_intron.push_back((*it1).Dig1);
                            else
                                fixed_intron.push_back((*it).Dig1);
                        }
                    }
                    else
                        cout<<"Error in plus_minus_flag"<<endl;

                }
                //修正完一行
                bed_line.block_bp_pair_number.clear();
                bed_line.DNA_block_begin_index.clear();
                int k=1;//k来索引new_intron
                bed_line.DNA_block_begin_index.push_back(0);
                bed_line.block_bp_pair_number.push_back(fixed_intron[0]-bed_line.DNA_comparison_begin_index);
                for(i=0;i<bed_line.comparison_block_number-2;i++){
                    bed_line.DNA_block_begin_index.push_back(fixed_intron[k]-bed_line.DNA_comparison_begin_index);
                    bed_line.block_bp_pair_number.push_back(fixed_intron[k+1]-fixed_intron[k]);
                    k+=2;
                }
                bed_line.DNA_block_begin_index.push_back(fixed_intron[k]-bed_line.DNA_comparison_begin_index);
                bed_line.block_bp_pair_number.push_back(bed_line.DNA_comparison_end_index-fixed_intron[k]);

                //把需要合并的合并

                vector<ll> tt1,tt2;//记录模块大小及开始位置
                tt1.clear();tt2.clear();
                tt1=bed_line.block_bp_pair_number;
                tt2=bed_line.DNA_block_begin_index;
                bed_line.block_bp_pair_number.clear();
                bed_line.DNA_block_begin_index.clear();
                int cc=0;//记录模块数
                bed_line.block_bp_pair_number.push_back(tt1[0]);
                bed_line.DNA_block_begin_index.push_back(tt2[0]);
                for(int i=0;i<bed_line.comparison_block_number;i++){
                    if(tt2[i]<=bed_line.DNA_block_begin_index[cc]+bed_line.block_bp_pair_number[cc]+maxIntronFix)
                        bed_line.block_bp_pair_number[cc]=tt2[i]+tt1[i]-bed_line.DNA_block_begin_index[cc];
                    else{
                        bed_line.block_bp_pair_number.push_back(tt1[i]);
                        bed_line.DNA_block_begin_index.push_back(tt2[i]);
                        cc++;
                    }
                }
                bed_line.comparison_block_number=cc+1;

               // Print_bed(bed_line,out);
               result_fixIntron.insert(bed_line);

            }
    }
    //输出排序、去重之后的文件s数据
    set<BEDDefine>::iterator fixed_it;
    for(fixed_it=result_fixIntron.begin(); fixed_it!=result_fixIntron.end(); fixed_it++)
    {
        Print_bed(*fixed_it,out);
    }

}

int main(int argc,char* argv[])
{
    int EvdFlag=0;//判断是否有证据文件，为1是有证据文件
    //******************设置短参数**********************************
    int oc;//短参数选项字符
    char* myopt;//短参数选项参数字符串
    while((oc=getopt(argc,argv,"a:b:c:d:"))!=EOF)
    {
        //短参数a跟输入文件路径，短参数b跟输出文件路径
        //短参数c跟证据文件路径，短参数d跟需要合并的最大intron值
        switch(oc)
        {
        case 'a':
            inputFilePath=optarg;//输入文件路径
            break;
        case 'b':
            //输出文件路径
            outputFilePath=optarg;
            break;
        case 'c':
            EvdFlag=1;
            myopt=optarg;
            evdFilePath.clear();
            evdFilePath.push_back(myopt);
            for(int i=0; i<argc-9; i++)
            {
                evdFilePath.push_back(argv[optind+i]);
            }
            break;
        case 'd':
            myopt=optarg;
            maxIntronFix=atoi(myopt);
            break;

        }
    }
    //*********************************************************

    all_index=0;
    freopen(inputFilePath,"r",stdin) ;

    Read() ;
    Solve() ;
    aa.clear();
    mode_one.clear();
    mode_two.clear();
    getmRNA_name.clear();
    sort(bed_data.begin(),bed_data.end());
    if(EvdFlag==1)
    {
        SecondStep();//有证据文件的
        //cout<<"end of secondstep()"<<endl;
        Evd_nag.clear();
        Evd_pos.clear();
        bed_data.clear();
        result_fixIntron.clear();
    }
    else
    {
        SecondStep_2();//没有证据文件的
        bed_data.clear();
        AccFu.clear();
        AccZheng.clear();
        DorFu.clear();
        DorZheng.clear();
        result_fixIntron.clear();
    }
    return 0;
}
