/* RzLi, Jan.21,2016
    This is a tool written at the 2nd year of my PhD research. I was switching from C style to C++ style. There might be some inconsistency in style(like dynamic arrays should be declared using "new *p").
 */

//Rongzhong Li, Jun.11,2011
//input pdb file should be in raw protein and rna pdb format. RNA_LENGTH and PROTEIN_LENGTH should be changed for larger molecules.
//the dihedral factor can be changed by searching "dihedral factor"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#define PI 3.1415926
#define RNA_LENGTH 2000
#define PROTEIN_LENGTH 2000	//for C_ALPHA
#define TIS_LENGTH RNA_LENGTH*3		//for P,B,S representation
#define CGM_LENGTH PROTEIN_LENGTH+TIS_LENGTH
#define PAIR_NUM TIS_LENGTH/2
#define HBOND 3.7
#define DIHEDRAL 30
#define COPLANARITY 8
#define TOLERANCE 1.3
#define TORSION "_torsions.dat"
#define N_T_HEAD "./gen_tis_torsions.pl -n -i "<<pdb<<".pdb -r "
#define P_T_HEAD "./gen_tis_torsions.pl -i "<<pdb<<".pdb -r "
#define STACKING "_stacks.dat"
#define S_HEAD "./gen_tis_stacking.pl -i "<<pdb<<".pdb "

using namespace std;

class coords
{
public:
    double x,y,z;
    coords()
    {x=y=z=0;}
};

class nucleotide
{
public:
    coords pho[4],sug[9],bas[20],C[10],O[10],N[10],norm,plane;
    char base[4];
    int pho_l,sug_l,bas_l,flag;
}rna[RNA_LENGTH];

class aminoacid
{
public:
    coords C_alpha;
    char name[4];
}protein[PROTEIN_LENGTH];

class pdbline
{
public:
    char line[2048],buffer[10];
    char atom[5],         atom_type[5],atom_name[4],chain_id[2],                             tail[20];
    int          atom_num,                                      resi_num;
    double                                                                x, y, z, a_occ,b_f;
    //PDB<| 4       7    1     4      1     3      1    1          4    4 8  8  8     6   6     19 |>
    //idx:      4       12           17           21          22        30 38 46 54    60  66
    pdbline()
    {
        atom_num=resi_num=x=y=z=a_occ=b_f=0;
    }
    int read(std::ifstream &raw_pdb)
    {
        raw_pdb.getline(line,2048);
        strncpy(atom,		line,	4);atom[4]		='\0';														//ATOM
        if(strcmp(atom,"ATOM"))return 1;
        if(!strcmp(atom,"END ")||!strcmp(atom,"TER "))return 2;
        if(!strcmp(atom,"ATOM"))
        {
            strncpy(buffer,		line+4,	7);	buffer[7]	='\0';atom_num	=atoi(buffer);								//atom_num
            strncpy(atom_type,	line+12,4);	atom_type[4]='\0';														//atom_type
            strncpy(atom_name,	line+17,3);	atom_name[3]='\0';														//atom_name
            strncpy(chain_id,	line+21,1);	chain_id[1]	='\0';														//chain_id
            strncpy(buffer,		line+22,4);	buffer[4]	='\0';resi_num	=atoi(buffer);								//resi_num
            strncpy(buffer,		line+30,8);	buffer[8]	='\0';x		=atof(buffer);									//x
            strncpy(buffer,		line+38,8);	buffer[8]	='\0';y		=atof(buffer);									//y
            strncpy(buffer,		line+46,8);	buffer[8]	='\0';z		=atof(buffer);									//z
            strncpy(buffer,		line+54,6);	buffer[6]	='\0';a_occ	=atof(buffer);									//atom_occupancy
            strncpy(buffer,		line+60,6);	buffer[6]	='\0';b_f	=atof(buffer);									//beta_factor
            strncpy(tail,			line+66,19);tail[19]	='\0';														//tail info
        }
        return 0;
    }
    void write()
    {
        cout		<<atom<<setw(7)<<atom_num<<setw(5)<<atom_type<<setw(4)<<atom_name<<setw(2)<<chain_id<<setw(4)<<resi_num<<setiosflags(ios::fixed)<<setprecision(3)
        <<setw(12)<<x<<setw(8)<<y<<setw(8)<<z<<setw(6)<<a_occ<<setw(6)<<b_f<<tail<<endl;
    }
    void fwrite(std::ofstream &out_pdb)
    {
        out_pdb	<<atom<<setw(7)<<atom_num<<setw(5)<<atom_type<<setw(4)<<atom_name<<setw(2)<<chain_id<<setw(4)<<resi_num<<setiosflags(ios::fixed)<<setprecision(3)
        <<setw(12)<<x<<setw(8)<<y<<setw(8)<<z//<<setw(6)<<a_occ<<setw(6)<<b_f<<tail
        <<endl;
    }
};

char pdbname[50];
int rna_length,protein_length,cgm_length,pairs,final_pairs;
int base_pair[PAIR_NUM][3];
int new_base_pair[PAIR_NUM][3];
int flags[TIS_LENGTH];
pdbline cgm[CGM_LENGTH];

//-----------------------Function Declarations-------------------------------------------
void readin_raw(char* raw);
void print_CGM();
void print_psf(char domain);
int base_pair_rule(int b1,int b2);
double distance(double x1,double y1,double z1,double x2,double y2,double z2);
double dist_ij(coords a,coords b);
double Hbond_dist_rule(int b1,int b2);
int find_base_pairing();
void normal_of_plane(int b);
double dihedral(int b1,int b2);
double coplanarity(int b1, int b2);
void check_pairing();
void print_base(int b);
void print_pair_info(int b1,int b2);
void gen_torsion_and_stacking(char *pdb);

//--------------------------------functions---------------------------------------

int main(int argn,char * argv[])
{
    if(argn!=2){cout <<"Format: "<<argv[0] <<" <pdb file>"<<endl;exit(-1);}
    int b1,b2,c=0;
    while(c<50){pdbname[c]=argv[1][c];c++;
        if(argv[1][c]=='.'&&argv[1][c+1]=='p'&&argv[1][c+2]=='d'&&argv[1][c+3]=='b')break;} //get the PDB name befor ".pdb"
    pdbname[c]='\0';
    ofstream out("./base_pairing.dat",ios::out);
    readin_raw(argv[1]);	//read in raw pdb file
    print_CGM();		//print out the CGM model
    print_psf('a');
    find_base_pairing();
    cout <<"Show pairing infomation of <n1><n2>.(Type '0' to exit)"<<endl;
    while (cin>>b1)
    {
        if (b1==0)break;
        cin >>b2;
        print_pair_info(b1-1,b2-1);
    }
    for(int i=0;i<final_pairs;i++)out <<new_base_pair[i][0]<<" "<<new_base_pair[i][1]<<endl;
    cout <<"Base pairing information is stored in './base_pairing.dat'."<<endl;
    gen_torsion_and_stacking(pdbname);
    cout <<"****************End of program****************"<<endl;
    out.clear();out.close();
}

void readin_raw(char* raw)	//read in raw pdb file
{
    ifstream raw_pdb(raw,ios::in);
    if(raw_pdb.fail()){cerr<<"Fail to open PDB file!"<<endl;exit(-1);}
    char line[2048],buffer[10],chain;
    char atom[5],         atom_type[5],atom_name[4],chain_id[2],                             tail[20];
    int          atom_num,                                      resi_num,                 			 resi_num_last=0;
    double                                                                x, y, z, a_occ,b_f;
    cerr <<"Reading raw PDB file...."<<endl;
    pdbline pdb;
    while(raw_pdb.good())				//reading raw pdb file
    {
        int token=pdb.read(raw_pdb);
        if(token==2)break;
        if(token==1)continue;
        strcpy(atom,pdb.atom);strcpy(atom_type,pdb.atom_type);strcpy(atom_name,pdb.atom_name);strcpy(chain_id,pdb.chain_id);strcpy(tail,pdb.tail);
        atom_num=pdb.atom_num;resi_num=pdb.resi_num;x=pdb.x;y=pdb.y;z=pdb.z;		//copy the info of pdbline for neater citation
        
        if(!strcmp(atom_type," CA "))
        {
            strcpy(protein[protein_length].name,atom_name);
            protein[protein_length].C_alpha.x=x;protein[protein_length].C_alpha.y=y;protein[protein_length].C_alpha.z=z;
            protein_length++;
            chain=chain_id[0];
            if(rna_length==1){resi_num_last=0;rna_length=0;}	//deal with the first residue which is mistaken as rna
            continue;}	//no need to read more if it's protein
        
        if(chain_id[0]!=chain&&resi_num_last<=resi_num)	//when the chain is not protein && forbids reading a second RNA chain
        {
            if(resi_num_last!=resi_num)					//do this only when meet a new residue
            {
                if(atom_name[0]=='A'||atom_name[0]=='U'||atom_name[0]=='C'||atom_name[0]=='G')
                {rna[rna_length].base[0]=atom_name[0];}
                else rna[rna_length].base[0]=atom_name[2];
                if			(rna[rna_length].base[0]=='A')strcpy(rna[rna_length].base,"ADE");
                else if	(rna[rna_length].base[0]=='U')strcpy(rna[rna_length].base,"URA");
                else if	(rna[rna_length].base[0]=='C')strcpy(rna[rna_length].base,"CYT");
                else if	(rna[rna_length].base[0]=='G')strcpy(rna[rna_length].base,"GUA");
                resi_num_last=resi_num;rna_length++;
            }
            pdb.resi_num=rna_length;
            //		pdb.write();//print the line with correct resi_num without skipping number
            //initialize rna array
            
            if(!strcmp(atom_type," OP1")||!strcmp(atom_type," OP2")||!strcmp(atom_type," O1P")||!strcmp(atom_type," O2P")||!strcmp(atom_type," P  ")||!strcmp(atom_type," OP3")||!strcmp(atom_type," O3P")
               )	//if it's phosphate
            {rna[rna_length-1].pho[rna[rna_length-1].pho_l].x=x;rna[rna_length-1].pho[rna[rna_length-1].pho_l].y=y;rna[rna_length-1].pho[rna[rna_length-1].pho_l].z=z;
                rna[rna_length-1].pho_l++;}								//P coords
            else if(atom_type[3]=='\''&&atom_type[0]!='H'&&atom_type[1]!='H')					//if it's heavy atom in suger
            {rna[rna_length-1].sug[rna[rna_length-1].sug_l].x=x;rna[rna_length-1].sug[rna[rna_length-1].sug_l].y=y;rna[rna_length-1].sug[rna[rna_length-1].sug_l].z=z;
                rna[rna_length-1].sug_l++;}								//S coords
            else if(atom_type[0]!='H'&&atom_type[1]!='H')										//heavy atom in base
            {rna[rna_length-1].bas[rna[rna_length-1].bas_l].x=x;rna[rna_length-1].bas[rna[rna_length-1].bas_l].y=y;rna[rna_length-1].bas[rna[rna_length-1].bas_l].z=z;
                rna[rna_length-1].bas_l++;								//B coords
                if(atom_type[1]=='C')
                {rna[rna_length-1].C[atom_type[2]-'0'].x=x;rna[rna_length-1].C[atom_type[2]-'0'].y=y;rna[rna_length-1].C[atom_type[2]-'0'].z=z;}
                if(atom_type[1]=='O')								//record important atom position of base plain
                {rna[rna_length-1].O[atom_type[2]-'0'].x=x;rna[rna_length-1].O[atom_type[2]-'0'].y=y;rna[rna_length-1].O[atom_type[2]-'0'].z=z;}
                if(atom_type[1]=='N')
                {rna[rna_length-1].N[atom_type[2]-'0'].x=x;rna[rna_length-1].N[atom_type[2]-'0'].y=y;rna[rna_length-1].N[atom_type[2]-'0'].z=z;}
            }
        }
    }
    
    if (protein_length)
        cerr <<protein_length<<" amino acids";
    if (rna_length)
    {
        if(protein_length)cerr <<" and ";
        cerr <<rna_length <<" nucleotides";}
    if(protein_length+rna_length==0)cerr<<"Nothing";
    cerr <<" imported."<<endl;
    raw_pdb.clear();raw_pdb.close();
    if(rna_length<10){cerr<<"There's no RNA data in the PDB file!"<<endl;exit(-1);}
    for(int l=0;l<rna_length;l++)
        normal_of_plane(l);
}

void print_CGM()	//print out the CGM model
{
    cout <<"****************Coordinates of CGM model****************"<<endl;
    double average_x,average_y,average_z;
    char outname[60];
    sprintf(outname,"./%s_CGM.pdb",pdbname);
    ofstream out(outname,ios::out);
    pdbline pdb;
    strcpy(pdb.atom,"ATOM");		//shared by all lines
    //-----------Protein-----------------------------------------------------------------
    if(protein_length)
        for(int l=0;l<protein_length;l++)
        {	strcpy(pdb.chain_id,"A");
            strcpy(pdb.atom_type," CA ");
            strcpy(pdb.atom_name,protein[l].name);
            pdb.resi_num=l+1;
            pdb.atom_num++;pdb.x=protein[l].C_alpha.x;pdb.y=protein[l].C_alpha.y;pdb.z=protein[l].C_alpha.z;
            cgm[cgm_length++]=pdb;
            pdb.fwrite(out);
        }
    //-----------RNA-----------------------------------------------------------------
    strcpy(pdb.chain_id,"B");
    for(int l=0;l<rna_length;l++)
    {int p,s,b;
        //cout <<l+1<<" "<<rna[l].base[0]<<endl;
        strcpy(pdb.atom_name,rna[l].base);pdb.resi_num=l+1;	//shared by this nucleotide;
        if(l)
        {
            for(p=0;p<rna[l].pho_l;p++){
                //cout <<rna[l].pho[p].x<<" "<<rna[l].pho[p].y<<" "<<rna[l].pho[p].z<<endl;
                average_x+=rna[l].pho[p].x;average_y+=rna[l].pho[p].y;average_z+=rna[l].pho[p].z;}
            strcpy(pdb.atom_type,"P ");
            pdb.atom_num++;pdb.x=average_x/p;pdb.y=average_y/p;pdb.z=average_z/p;
            cgm[cgm_length++]=pdb;
            pdb.fwrite(out);
            //out <<setw(8)<<average_x/p<<setw(8)<<average_y/p<<setw(8)<<average_z/p<<endl;
            average_x=average_y=average_z=0;
        }
        for(s=0;s<rna[l].sug_l;s++){
            //cout <<rna[l].sug[s].x<<" "<<rna[l].sug[s].y<<" "<<rna[l].sug[s].z<<endl;
            average_x+=rna[l].sug[s].x;average_y+=rna[l].sug[s].y;average_z+=rna[l].sug[s].z;}
        strcpy(pdb.atom_type,"S ");
        pdb.atom_num++;pdb.x=average_x/s;pdb.y=average_y/s;pdb.z=average_z/s;
        cgm[cgm_length++]=pdb;
        pdb.fwrite(out);
        //out <<setw(8)<<average_x/s<<setw(8)<<average_y/s<<setw(8)<<average_z/s<<endl;
        average_x=average_y=average_z=0;
        for(b=0;b<rna[l].bas_l;b++){
            //cout <<rna[l].bas[b].x<<" "<<rna[l].bas[b].y<<" "<<rna[l].bas[b].z<<endl;
            average_x+=rna[l].bas[b].x;average_y+=rna[l].bas[b].y;average_z+=rna[l].bas[b].z;}
        strcpy(pdb.atom_type,"B ");
        pdb.atom_num++;pdb.x=average_x/b;pdb.y=average_y/b;pdb.z=average_z/b;
        cgm[cgm_length++]=pdb;
        pdb.fwrite(out);
        //out <<setw(8)<<average_x/b<<setw(8)<<average_y/b<<setw(8)<<average_z/b<<endl;
        average_x=average_y=average_z=0;
    }
    cout <<"Coarse-Grained Model(CGM) coordinates are stored in './"<<pdbname<<"_CGM.pdb'."<<endl;
    out.clear();out.close();
    cout <<"***********************Sequence**************************"<<endl;
    cout <<rna_length<<" nucleotides:"<<endl;
    for(int l=0;l<rna_length;l++)cout <<rna[l].base[0];
    cout <<endl;
}

void print_psf(char domain)	//print out the PSF info
{
    char outname[60];int bond=0,bonds[CGM_LENGTH][2];
    int cgm_l,start,end,shift;
    sprintf(outname,"./%s_CGM.psf",pdbname);
    ofstream out(outname,ios::out);
    out	<<"PSF CMAP"<<endl<<endl
    <<setw(8)<<0<<" !NTITLE"<<endl<<endl;
    if(domain=='r'){start=protein_length;end=cgm_length;}
    else if(domain=='p'){start=0;end=protein_length;}
    else {start=0;end=cgm_length;}
    //--------------atoms----------------------------------
    out	<<setw(8)<<end-start<<"!NATOM"<<endl;
    for(int cgm_l=start;cgm_l<end;cgm_l++)
        out <<setw(8)<<cgm_l+1-start<<setw(5)<<cgm[cgm_l].chain_id<<' '<<setiosflags(ios::left)<<setw(4)<<cgm[cgm_l].resi_num
        <<resetiosflags(ios::left)<<setiosflags(ios::right)<<setw(4)<<cgm[cgm_l].atom_name<<setw(5)<<cgm[cgm_l].atom_type<<setw(5)<<cgm[cgm_l].atom_type
        <<setiosflags(ios::fixed)<<setprecision(6)<<setw(12)<<0.000000<<setprecision(4)<<setw(14)<<12.0110<<setw(12)<<0<<endl<<resetiosflags(ios::right);
    out <<endl;
    //---------------bonds-------------------------------
    if(domain!='r')//p or a
        for(int l=0;l<protein_length-1;l++)
        {bonds[bond][0]=l+1;bonds[bond][1]=l+2;bond++;}
    if(domain!='p')//r or a
    {if(domain=='r')shift=0;else shift=protein_length;
        for(int l=0;l<rna_length;l++)
        {
            if(l)
            {bonds[bond][0]=l*3-2+shift;bonds[bond][1]=l*3+shift;bond++;
                bonds[bond][0]=l*3+shift;bonds[bond][1]=l*3+1+shift;bond++;}
            bonds[bond][0]=l*3+1+shift;bonds[bond][1]=l*3+2+shift;bond++;
        }
    }
    out	<<setw(8)<<bond<<" !NBOND: bonds"<<endl;
    for(int b=0;b<bond;b++)
    {if(!(b%4)&&b)out<<endl;
        out <<setw(8)<<bonds[b][0]<<setw(8)<<bonds[b][1];
    }
    out.clear();out.close();
}


int base_pair_rule(int b1,int b2)//compare with know base pairing rule. 0 for fail,1 for normal,2 for others.
{
    char b1_base=rna[b1].base[0],b2_base=rna[b2].base[0];
    switch(b1_base){
        case 'A':	if(b2_base=='U')return 1;else if(b2_base=='C')return 2;else if(b2_base=='G')return 3;else return 0;
        case 'G':	if(b2_base=='C')return 1;else if(b2_base=='U')return 2;else if(b2_base=='A')return 3;else return 0;
        case 'C':	if(b2_base=='G')return 1;else if(b2_base=='A')return 2;else return 0;
        case 'U':	if(b2_base=='A')return 1;else if(b2_base=='G')return 2;else return 0;
        default:	return 0;
    }}


double Hbond_dist_rule(int b1,int b2)		//calculate the exact length of H-bond, suppose it exists.
{
    if (b1==b2)return 0;
    int swap;
    char base1=rna[b1].base[0],base2=rna[b2].base[0];//check the base name
    //---------------below is a list of known base pairing rules--------------------
    if(base1=='A'&&base2=='U'||base2=='A'&&base1=='U')
    {if(base2=='A'&&base1=='U'){swap=b1;b1=b2;b2=swap;}
        return min(dist_ij(rna[b1].N[6],rna[b2].O[4]),
                   min(dist_ij(rna[b1].N[6],rna[b2].O[2]),
                       min(dist_ij(rna[b1].N[1],rna[b2].O[3]),
                           dist_ij(rna[b1].N[7],rna[b2].O[3])
                           )
                       )
                   );
    }
    else if(base1=='G'&&base2=='U'||base2=='G'&&base1=='U')
    {if(base2=='G'&&base1=='U'){swap=b1;b1=b2;b2=swap;}
        return	min(dist_ij(rna[b1].O[6],rna[b2].N[3]),
                    min(dist_ij(rna[b1].N[1],rna[b2].O[2]),
                        dist_ij(rna[b1].N[1],rna[b2].O[4])
                        )
                    );
    }
    else if(base1=='G'&&base2=='C'||base2=='G'&&base1=='C')
    {if(base2=='G'&&base1=='C'){swap=b1;b1=b2;b2=swap;}
        return	min(dist_ij(rna[b1].O[6],rna[b2].N[4]),
                    min(dist_ij(rna[b1].N[1],rna[b2].N[3]),
                        min(dist_ij(rna[b1].N[2],rna[b2].O[2]),
                            min(dist_ij(rna[b1].N[1],rna[b2].O[2]),
                                dist_ij(rna[b1].N[2],rna[b2].N[3])
                                )
                            )
                        )
                    );
    }
    else if(base1=='A'&&base2=='C'||base2=='A'&&base1=='C')
    {if(base2=='A'&&base1=='C'){swap=b1;b1=b2;b2=swap;}
        return min(dist_ij(rna[b1].N[6],rna[b2].N[3]),
                   min(dist_ij(rna[b1].N[7],rna[b2].N[4]),
                       dist_ij(rna[b1].N[1],rna[b2].N[4])
                       )
                   );
    }
    else return HBOND*10;//if no such rule, return a large number
}


double dist_ij(coords a,coords b)
{
    return distance(a.x,a.y,a.z,b.x,b.y,b.z);
}

double distance(double x1,double y1,double z1,double x2,double y2,double z2)
{return sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) );}



void normal_of_plane(int b)//calculate the normal vector of base plane
{
    double x1,y1,z1,x2,y2,z2,r;
    x1=rna[b].C[2].x-rna[b].C[4].x;x2=rna[b].C[6].x-rna[b].C[4].x;
    y1=rna[b].C[2].y-rna[b].C[4].y;y2=rna[b].C[6].y-rna[b].C[4].y;
    z1=rna[b].C[2].z-rna[b].C[4].z;z2=rna[b].C[6].z-rna[b].C[4].z;
    rna[b].norm.x=-y2*z1+y1*z2;
    rna[b].norm.y= x2*z1-x1*z2;
    rna[b].norm.z=-x2*y1+x1*y2;
    r=sqrt(rna[b].norm.x*rna[b].norm.x+rna[b].norm.y*rna[b].norm.y+rna[b].norm.z*rna[b].norm.z);
    rna[b].norm.x/=r;
    rna[b].norm.y/=r;
    rna[b].norm.z/=r;
}

double coplanarity(int b1, int b2)
{
    double cosine,x1,y1,z1,x2,y2,z2;
    if (b1==b2)return 0;
    x1=rna[b1].C[4].x-rna[b2].C[4].x;
    y1=rna[b1].C[4].y-rna[b2].C[4].y;
    z1=rna[b1].C[4].z-rna[b2].C[4].z;
    if(dihedral(b1,b2)<PI/2)
    {x2=rna[b1].norm.x+rna[b2].norm.x;y2=rna[b1].norm.y+rna[b2].norm.y;z2=rna[b1].norm.z+rna[b2].norm.z;}
    else
    {x2=rna[b1].norm.x-rna[b2].norm.x;y2=rna[b1].norm.y-rna[b2].norm.y;z2=rna[b1].norm.z-rna[b2].norm.z;}
    cosine=(x1*x2+y1*y2+z1*z2)/sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2));
    return fabs(cosine);
}

double dihedral(int b1,int b2)//base pairing should be in a plane
{
    double cosine;
    cosine=(rna[b1].norm.x*rna[b2].norm.x+rna[b1].norm.y*rna[b2].norm.y+rna[b1].norm.z*rna[b2].norm.z)/sqrt((rna[b1].norm.x*rna[b1].norm.x+rna[b1].norm.y*rna[b1].norm.y+rna[b1].norm.z*rna[b1].norm.z) *(rna[b2].norm.x*rna[b2].norm.x+rna[b2].norm.y*rna[b2].norm.y+rna[b2].norm.z*rna[b2].norm.z));
    return acos(cosine);
}

int find_base_pairing()	//find base pairing with judgement of distance, pairing rules and angle between base plane
//later it calls a checking function to get a refined list
{
    double tolerance;
    for(int p=0;p<PAIR_NUM;p++){base_pair[p][0]=0;base_pair[p][1]=0;base_pair[p][2]=0;}
    cerr <<"****************List of base pairing****************"<<endl;
    for(int b1=0;b1<rna_length;b1++)
        for(int b2=b1+2;b2<rna_length;b2++)
        {
            if(base_pair[pairs-1][0]==b1&&base_pair[pairs-1][1]==b2+2)tolerance=TOLERANCE;	//looser restriction
            if	(
                 base_pair_rule(b1,b2)		//pairing rule
                 &&Hbond_dist_rule(b1,b2)<=HBOND*tolerance	//H-bond distance
                 &&(dihedral(b1,b2)/PI*180<DIHEDRAL*tolerance||dihedral(b1,b2)/PI*180>(180-DIHEDRAL*tolerance))//dihedral between base planes
                 &&coplanarity(b1,b2)/PI*180<COPLANARITY*tolerance//coplanarity of base planes
                 )
            {print_pair_info(b1,b2);
                base_pair[pairs][0]=b1+1;base_pair[pairs][1]=b2+1;base_pair[pairs][2]=base_pair_rule(b1,b2);
                pairs++;
            }
            tolerance=1;
        }
    cout <<pairs<<" base pairing found!"<<endl;
    check_pairing();
}

void check_pairing()//delete the pairing which shares bases
{
    cout <<"****************checking*********************"<<endl;
    int del=0;
    //int newlist=0;
    //new_base_pair[0][0]=base_pair[0][0];new_base_pair[0][1]=base_pair[0][1];newlist++;//copy right pairing to a new list
    for(int p=1;p<pairs;p++)
    {//new_base_pair[newlist][0]=base_pair[p][0];new_base_pair[newlist][1]=base_pair[p][1];newlist++;
        if(!base_pair[p][2])continue;
        int inc_1l,inc_1r,inc_2l,inc_2r;
        inc_1l=base_pair[p][0]-base_pair[p-1][0];inc_1r=base_pair[p+1][0]-base_pair[p][0];
        inc_2l=base_pair[p][1]-base_pair[p-1][1];inc_2r=base_pair[p+1][1]-base_pair[p][1];
        int b1=base_pair[p][0]-1,b2=base_pair[p][1]-1;
        
        
        if((abs(inc_1l)>1&&abs(inc_1r)>1)||(abs(inc_2l)>1&&abs(inc_2r)>1))
        {
            cout <<"*******Delete ";print_base(b1);print_base(b2);cout <<"because it's ISOLATED.*******"<<endl;
            base_pair[p][0]=base_pair[p-1][0];base_pair[p][1]=base_pair[p-1][1];base_pair[p][2]=0;del++;
        }
        
    }
    ofstream out("./pairing0.dat",ios::out);
    final_pairs=0;
    cout <<"****************Final List*****************"<<endl;
    for(int p=0;p<pairs+1;p++,final_pairs++)
    {
        if(!base_pair[p][2]){final_pairs--;continue;}
        int b1=base_pair[p][0]-1,b2=base_pair[p][1]-1;
        print_base(b1);print_base(b2);
        for(int j=0;j<3;j++)new_base_pair[final_pairs][j]=base_pair[p][j];//store new list
        out <<new_base_pair[final_pairs][0]<<"\t"<<new_base_pair[final_pairs][1]<<endl;
        if(base_pair_rule(b1,b2)==2)cout <<"*";
        else if(base_pair_rule(b1,b2)==3)cout <<"**";cout <<endl;
    }
    cout <<final_pairs<<" base pairing found!"<<endl;
    out.clear();out.close();
}

void print_base(int b)
{
    cout <<b+1<<rna[b].base[0]<<"\t";
}

void print_pair_info(int b1,int b2)
{
    print_base(b1);print_base(b2);cout <<"Type: "<<base_pair_rule(b1,b2)<<"\tBond: "<<Hbond_dist_rule(b1,b2)<<"\tDihedral: "<<dihedral(b1,b2)*180/PI<<"\tCoplanarity: "<<coplanarity(b1,b2)/PI*180<<endl;
}


//----------------------------------------------------------
void gen_torsion_and_stacking(char *pdb)
{
    ofstream setupfile("./setup_data.sh",ios::out);
    setupfile	<<"#!/bin/bash"<<endl<<endl
    <<"./gen_tis_input.pl -i "	<<pdb	<<".pdb -b >\t../DATA/"					<<pdb<<"_bonds.dat"	<<endl
    <<"./gen_tis_input.pl -i "	<<pdb	<<".pdb -a >\t../DATA/"					<<pdb<<"_angles.dat"<<endl
    <<"./gen_tis_input.pl -i "	<<pdb	<<".pdb -v >\t../DATA/"					<<pdb<<"_nb.dat"	<<endl
    <<"./gen_tis_input.pl -i "	<<pdb	<<".pdb -e >\t../DATA/"					<<pdb<<"_elec.dat"	<<endl
    <<"echo \"nbead "			<<protein_length+rna_length*3-1<<"\" > ../DATA/"<<pdb<<"_init.xyz"	<<endl
    <<"./gen_tis_input.pl -i "	<<pdb	<<".pdb -x >>\t../DATA/"				<<pdb<<"_init.xyz"	<<endl
    <<"./gen_tis_input.pl -i "	<<pdb	<<".pdb -x >\t../DATA/"					<<pdb<<"_nat.xyz"	<<endl<<endl;
    string factor[3]={" -f 0.5 "," "," -f 1 "};//change dihedral factor here: rna stem, loop, protein
    for(int l=0;l<rna_length;l++)flags[l]=0;
    for(int j=0;j<final_pairs;j++)
    {flags[new_base_pair[j][0]-1]=1;flags[new_base_pair[j][1]-1]=1;}//mark flags
    //torsion part
    setupfile <<"echo \"ntorsion "<<protein_length-3+rna_length*2-4<<"\" >\t../DATA/"<<pdb<<TORSION<<endl;
    //-----protein-----
    for(int i=0;i<protein_length-2;i++)
    {
        if(i==0){setupfile <<P_T_HEAD<<" "<<i+1;continue;}
        setupfile<<":"<<i+1;
        if(i==protein_length-3)
            setupfile <<factor[2]<<">>\t../DATA/"<<pdb<<TORSION<<endl;
    }
    //-----rna----------
    for(int i=0;i<rna_length-1;i++)
    {
        if(i==0||flags[i-1]!=flags[i])setupfile <<N_T_HEAD<<" "<<i+1+protein_length;
        else setupfile<<":"<<i+1+protein_length;
        if(i==rna_length-2||flags[i+1]!=flags[i])
        {if(flags[i])	setupfile <<factor[flags[i]];
        else 			setupfile <<factor[flags[i]];
            setupfile <<">>\t../DATA/"<<pdb<<TORSION<<endl;}
    }
    setupfile <<endl;
    //stacking part
    int nstack=0,inc_l,inc_r;
    for(int j=0;j<final_pairs-1;j++)
    {
        inc_l=new_base_pair[j+1][0]-new_base_pair[j][0];
        inc_r=new_base_pair[j+1][1]-new_base_pair[j][1];
        if(abs(inc_l)==1&&abs(inc_r)==1)
        {
            nstack++;
        }
    }
    setupfile <<"echo \"nstack "<<nstack<<"\" >\t../DATA/"<<pdb<<STACKING<<endl;
    for(int j=0;j<final_pairs-1;j++)
    {
        inc_l=new_base_pair[j+1][0]-new_base_pair[j][0];
        inc_r=new_base_pair[j+1][1]-new_base_pair[j][1];
        if(abs(inc_l)==1&&abs(inc_r)==1)
        {
            setupfile	<<S_HEAD
            <<" -a "<<new_base_pair[j][0]+protein_length<<" -b "<<new_base_pair[j+1][0]+protein_length<<" -c "<<new_base_pair[j+1][1]+protein_length<<" -d "<<new_base_pair[j][1]+protein_length
            <<" >>\t../DATA/"<<pdb<<STACKING<<endl;
        }
    }
    setupfile.clear();setupfile.close();
    system("chmod u+x setup_data.sh");
    cout<<"Setup information is stored in './setup_data.sh'."<<endl;
}
