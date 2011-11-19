#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <functional>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>

#ifdef INDEP_PROGRAM
#define PRINTFUNCTION printf
#else
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#define PRINTFUNCTION Rprintf
#endif

using namespace std;

int scores[576] = {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4,
-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4,
-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4,
-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4,
0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,
-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4,
-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4,
0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4,
-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4,
-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4,
-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4,
-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4,
-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4,
-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4,
-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4,
1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4,
0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4,
-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4,
-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4,
0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4,
-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4,
-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4,
0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4,
-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1};
char acids[24] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','-'};

map<string, int> pair2score;
void print_blosum62();
void random_check_blosum62();
int get_score_seqs(string & str1, string & str2);
vector<string> string_tokenize(const string& str, const string& delimiters = " \t\n\r", bool skip_empty = true);

vector<pair<string, double> > score_matrix;
map<string, int> geneName2rank;

map<string, int> predefined;
map<string, int> predefined_seq;
map<string, string> id2seq;
vector<string> total_poteins;

double getES(vector<int>& indexes);
double getES_curves(vector<int>& indexes);
double getES();
double getES_curves();
double test_a_protein1(string& protein_name, double & ES);
double test_a_protein2(string& protein_name, double & ES);
int p_times=1000;
//ofstream Log("log.txt");
bool print_score_matrix();
int print_curves = 0;
ofstream CURVE;
string line1;
string line2;
string int2str(int value);
string double2str(double value);
void ToUpperString(string &str);
void aseb_sites(string input1, string input2, string input3, string output, int Permuationtimes);
void aseb_protein(string input1, string input2, string input3, string output, int Permuationtimes);
#ifdef INDEP_PROGRAM
int main(int argc, char *argv[])
{
  if(argc != 7){
     cerr<<"aseb sites <background_sites> <predefined_sites> <sites_to_test> <output> <permutation_times>"<<endl;
     cerr<<"aseb proteins <background_sites> <predefined_sites> <proteins_to_test> <output> <permutation_times>"<<endl;
     exit(1);
  }
  ///////////////////
  //print_curves = 1;
  //CURVE.open("test.txt");
  ///////////////////
  string input1 = argv[2];
  string input2 = argv[3];
  string input3 = argv[4];
  string output = argv[5];
  int Permuationtimes;
  sscanf(argv[6], "%d", &Permuationtimes);
  if(string(argv[1]) == "sites"){
     aseb_sites(input1, input2, input3, output, Permuationtimes);
  }
  if(string(argv[1]) == "proteins"){
     aseb_protein(input1, input2, input3, output, Permuationtimes);
  }
}
#else
#endif
void aseb_sites(string input1, string input2, string input3, string output, int Permuationtimes){
  int permutation_kind = 0;
  p_times = Permuationtimes;
  for(int i=0; i < 25; i++){
      for(int j=0; j < 25; j++){
          char tmp[10];
          sprintf(tmp, "%c_%c", acids[i], acids[j]);
          pair2score[tmp] = scores[i*24+j];
      }
  }
  ifstream in1(input1.c_str());
  if(!in1){
     PRINTFUNCTION("Can not open %s\n", input1.c_str());
     return;
  }
  ifstream in2(input2.c_str());
  if(!in2){
     PRINTFUNCTION("Can not open %s\n", input2.c_str());
     return;
  }
  ifstream in3(input3.c_str());
  if(!in3){
     PRINTFUNCTION("Can not open %s\n", input3.c_str());
     return;
  }
  ofstream out(output.c_str());
  if(!out){
     PRINTFUNCTION("Can not open %s\n", output.c_str());
     return;
  }

  char buffer[10000 + 1];
  
  string id = "";
  string seq = "";
  while(!in1.eof()){
      in1.getline(buffer, 1000000);
      string tmp = buffer;
      if(tmp[tmp.size() - 1] == '\r'){
         tmp[tmp.size() - 1] = '\0';
      }
      if(tmp[0] == '>'){
      	 vector<string> tokens = string_tokenize(tmp.substr(1));
      	 if(id != ""){
             total_poteins.push_back(id);
      	     id2seq[id] = seq;
      	 }
      	 id = tokens[0];
         seq = "";
      }else{
      	 ToUpperString(tmp);
         seq += tmp;
      }
  }
  if(id != ""){
     total_poteins.push_back(id);
     id2seq[id] = seq;
  }
  id = "";
  seq = "";

  while(!in2.eof()){
      in2.getline(buffer, 1000000);
      string tmp = buffer;
      if(tmp[tmp.size() - 1] == '\r'){
         tmp[tmp.size() - 1] = '\0';
      }
      if(tmp[0] == '>'){
      	 vector<string> tokens = string_tokenize(tmp.substr(1));
      	 if(id != ""){
      	 	  predefined[id] = 1;
            predefined_seq[seq] = 1;
            id2seq[id] = seq;
      	 }
      	 id = tokens[0];
         seq = "";
      }else{
      	 ToUpperString(tmp);
         seq += tmp;
      }
  }
  if(id != ""){
     predefined[id] = 1;
     predefined_seq[seq] = 1;
     id2seq[id] = seq;
  }
  id = "";
  seq = "";
  
  int count = 0;

  while(!in3.eof()){
      in3.getline(buffer, 1000000);
      string tmp = buffer;
      if(tmp[tmp.size() - 1] == '\r'){
         tmp[tmp.size() - 1] = '\0';
      }
      if(tmp[0] == '>'){
      	 vector<string> tokens = string_tokenize(tmp.substr(1));
      	 if(id != ""){
            id2seq[id] = seq;
            double p_value = 1;
            double ES_value = 0;
            if(permutation_kind == 0){
               p_value = test_a_protein1(id, ES_value);
            }else{
               p_value = test_a_protein2(id, ES_value);
            }
            out<<id<<"\t"<<ES_value<<"\t"<<p_value<<"\n";
            count++;
            if(count%100==0){
               PRINTFUNCTION("\rprocessed %d sites", count);
               #ifdef INDEP_PROGRAM
               #else
                 R_FlushConsole();
               #endif
            }
      	 }
      	 id = tokens[0];
         seq = "";
      }else{
      	 ToUpperString(tmp);
         seq += tmp;
      }
  }
  if(id != ""){
     id2seq[id] = seq;
     double p_value = 1;
     double ES_value = 0;
     if(permutation_kind == 0){
        p_value = test_a_protein1(id, ES_value);
     }else{
        p_value = test_a_protein2(id, ES_value);
     }
     out<<id<<"\t"<<ES_value<<"\t"<<p_value<<"\n";
     count++;
     if(count%100==0){
        PRINTFUNCTION("\rprocessed %d sites", count);
        #ifdef INDEP_PROGRAM
        #else
        R_FlushConsole();
        #endif
     }
  }
  id = "";
  seq = "";
  ///////////////////////////////////////////////////
  PRINTFUNCTION("\rprocessed %d sites\n", count);
  #ifdef INDEP_PROGRAM
  #else
  R_FlushConsole();
  #endif
  CURVE.close();
  out.close();
}
void aseb_protein(string input1, string input2, string input3, string output, int Permuationtimes){
  ofstream tmpOut((output+".tmp").c_str());
  if(!tmpOut){
     //cerr<<"Can not open "<<(output+".tmp")<<endl;
     PRINTFUNCTION("Can not open %s \n", (output+".tmp").c_str());
     return;
  }
  ifstream in3(input3.c_str());
  if(!in3){
     //cerr<<"Can not open "<<input3<<endl;
     PRINTFUNCTION("Can not open %s \n", (input3).c_str());
     return;
  }
  char buffer[1000000 + 1];
  string id = "";
  string seq = "";
  while(!in3.eof()){
      in3.getline(buffer, 1000000);
      string tmp = buffer;
      if(tmp[tmp.size() - 1] == '\r'){
         tmp[tmp.size() - 1] = '\0';
      }
      if(tmp[0] == '>'){
      	 vector<string> tokens = string_tokenize(tmp.substr(1));
      	 if(id != ""){
      	 	  string sequence = "----"+seq+"----";
            for(int i=8; i < (int)sequence.size()-8;i++){
                if(sequence[i] == 'K'){
                   string site_seq = sequence.substr(i-8, 17);
                   tmpOut<<">"<<id<<"_"<<i-4+1<<"\n"<<site_seq<<"\n";
                }
            }
      	 }
      	 id = tokens[0];
         seq = "";
      }else{
      	 ToUpperString(tmp);
         seq += tmp;
      }
  }
  if(id != ""){
     string sequence = "----"+seq+"----";
     for(int i=8; i < (int)sequence.size()-8;i++){
         if(sequence[i] == 'K'){
            string site_seq = sequence.substr(i-8, 17);
            tmpOut<<">"<<id<<"_"<<i-4+1<<"\n"<<site_seq<<"\n";
         }
     }
  }
  tmpOut.close();
  aseb_sites(input1, input2, (output+".tmp"), output, Permuationtimes);
}
inline bool bigThan(const pair<string, double> &a, const pair<string, double> &b){
   if(a.second - b.second > 0.001){
      return true;
   }
   return false;
}
bool generate_score_matrix(string& protein_name){
   score_matrix.clear();
   for(int i=0; i < (int)total_poteins.size(); i++){
       int score = get_score_seqs(id2seq[protein_name], id2seq[total_poteins[i]]);
       score_matrix.push_back(pair<string, int>(total_poteins[i], score));
   }
   //make_heap(score_matrix.begin(), score_matrix.end(), bigThan);
   //score_matrix.begin(), score_matrix.end(), bigThan);
   //sort(score_matrix.begin(), score_matrix.end(), bigThan);
   stable_sort(score_matrix.begin(), score_matrix.end(), bigThan);
   geneName2rank.clear();
   double max_pos = 1;
   double min_neg = -1;
   if(score_matrix[0].second > 0){
      max_pos = score_matrix[0].second;
   }
   if(score_matrix[score_matrix.size()-1].second < 0){
      min_neg = score_matrix[score_matrix.size()-1].second;
   }
   for(int i=0; i < (int)score_matrix.size(); i++){
       if(score_matrix[i].second > 0){
          score_matrix[i].second = score_matrix[i].second/max_pos;
       }
       if(score_matrix[i].second < 0){
          score_matrix[i].second = score_matrix[i].second/min_neg;
       }
       geneName2rank[score_matrix[i].first] = i;
   }
   return true;
}
bool print_score_matrix(){
   for(int i=0; i < (int)score_matrix.size(); i++){
       //Log<<score_matrix[i].first<<": "<<score_matrix[i].second<<endl;
   }
   //Log<<"---------------------------------"<<endl;
   return true;
}
void get_random_indexes(vector<int>& indexes){
   map<int, int> wheter_have;
   int max = score_matrix.size();
   int size = predefined.size();
   indexes.clear();
   indexes.reserve(size);
   int rand_upper = ((RAND_MAX-max+1)/max)*max-1+max;
   while((int)indexes.size() < size){
     int rand_num1 = rand();
     int rand_num2 = rand_num1%max;
     while((rand_num1 > rand_upper)||(wheter_have.count(rand_num2) != 0)){
         rand_num1 = rand();
         rand_num2 = rand_num1%max;
     }
     indexes.push_back(rand_num2);
     wheter_have[rand_num2]=1;
   }
   //sort(indexes.begin(), indexes.end());
   /*Log<<"--------indexes"<<endl;
   for(int i=0; i < (int)indexes.size(); i++){
       Log<<indexes[i]<<endl;
   }
   Log<<endl;*/
}
void get_random_indexes2(vector<int>& random_index){
   random_index.clear();
   map<int, int> exist;
   int max = total_poteins.size();
   int rand_upper = ((RAND_MAX-max+1)/max)*max-1+max;
   while(random_index.size() < predefined.size()){
      int rand_num1 = rand();
      int rand_num2 = rand_num1%max;
      while((rand_num1 > rand_upper)||(exist.count(rand_num2) != 0)){
          rand_num1 = rand();
          rand_num2 = rand_num1%max;
      }
      exist[rand_num2] = 1;
      random_index.push_back(geneName2rank[total_poteins[rand_num2]]);
      //random_index.push_back(geneName2rank[score_matrix[rand_num2].first]);
   }
   sort(random_index.begin(), random_index.end());
}
void swap(char *pm,char *pn)
{
   char temp;
   temp=*pm;
   *pm=*pn;
   *pn=temp;
}
string get_random_str(string &str1){
   string rand_str = str1;
   int size = str1.size();
   vector<int> indexes;
   for(int i=0; i < size/2; i++){
       indexes.push_back(i);
   }
   for(int i = size/2+1; i < size; i++){
       indexes.push_back(i);
   }
   random_shuffle(indexes.begin(), indexes.end());
   for(int i=0; i < size/2; i++){
       rand_str[i] = str1[indexes[i]];
       //cerr<<i<<" "<<indexes[i]<<endl;
   }
   for(int i=size/2+1; i < size; i++){
       rand_str[i] = str1[indexes[i-1]];
       //cerr<<i<<" "<<indexes[i-1]<<endl;
   }
   return rand_str;
}
double test_a_protein1(string& protein_name, double & ES){
   srand(12345);
   if(predefined_seq.count(id2seq[protein_name]) != 0){
      return ((double)1/p_times);
   }
   generate_score_matrix(protein_name);
   print_score_matrix();
   double ES_protein;
   if(print_curves > 0){
      line1 = protein_name + "\t";
      line2 = protein_name + "\t";
      ES_protein = getES_curves();
   }else{
      ES_protein = getES();
   }
   ES = ES_protein;
   double ES_random;
   int p_count=1;
   for(int i=1; i < p_times; i++){
       vector<int> random_index;
       get_random_indexes(random_index);
       sort(random_index.begin(), random_index.end());
       /*for(int j=0; j < (int)random_index.size(); j++){
           Log<<random_index[j]<<endl;
       }*/
       ES_random = getES(random_index);
       if(ES_random > ES_protein+0.0000001){
          p_count++;
       }
   }
   if(print_curves > 0){
       line1 += double2str(ES_protein)+"\t"+double2str((double)p_count/p_times)+"\n";
       line2 += double2str(ES_protein)+"\t"+double2str((double)p_count/p_times)+"\n";
       CURVE<<line1;
       CURVE<<line2;
   }
   //Log<<"ES 0:"<<ES_protein<<endl;
   return ((double)p_count/p_times);
}
double test_a_protein2(string& protein_name, double &ES){
   srand(12345678);
   if(predefined_seq.count(id2seq[protein_name]) != 0){
      return ((double)1/p_times);
   }
   generate_score_matrix(protein_name);
   double ES_protein;
   if(print_curves > 0){
      line1 = protein_name + "\t";
      line2 = protein_name + "\t";
      ES_protein = getES_curves();
   }else{
      ES_protein = getES();
   }
   ES = ES_protein;
   double ES_random;
   int p_count=1;
   string sequence = id2seq[protein_name];
   for(int i=1; i < p_times; i++){
       id2seq[protein_name] = get_random_str(sequence);
       generate_score_matrix(protein_name);
       id2seq[protein_name] = sequence;
       ES_random = getES();
       if(ES_random > ES_protein+0.00001){
          p_count++;
       }
   }
   //Log<<"ES 0:"<<ES_protein<<endl;
   return ((double)p_count/p_times);
}
void random_check_blosum62(){
    int test_time = 10;
    while(test_time > 0){
         int i = rand()%24;
         int j = rand()%24;
         char tmp[10];
         sprintf(tmp, "%c_%c", acids[i], acids[j]);
         cerr<<acids[i]<<"_"<<acids[j]<<": "<<pair2score[tmp]<<endl;
         test_time--;
    }
}
void print_blosum62(){
    for(int i=0; i < 24; i++){
        cerr<<acids[i]<<" ";
    }
    cerr<<endl;
    for(int i=0; i < 576; i++){
        cerr<<scores[i]<<" ";
        if((i+1)%24==0){
           cerr<<acids[((i+1)-(i+1)%24)/24-1]<<"\n";
        }
    }
    cerr<<endl;
}
int get_score_seqs(string & str1, string & str2){
    int score = 0;
    if(str1.size() != str2.size()){
       //cerr<<str1<<endl;
       //cerr<<str2<<endl;
       //cerr<<"Different length"<<endl;
       PRINTFUNCTION("Different length!\n");
       return -1;
    }
    for(int i=0; i < (int)str1.size(); i++){
        char tmp[10];
        char chr1 = str1[i];
        char chr2 = str2[i];
        if(chr1=='U'){
           chr1='-';
        }
        if(chr2=='U'){
           chr2='-';
        }
        sprintf(tmp, "%c_%c", chr1, chr2);
        if(pair2score.count(tmp)==0){
           //cerr<<tmp<<endl;
           //cerr<<str1<<endl;
           //cerr<<str2<<endl;
           //cerr<<"Contains unrecognizable character"<<endl;
           PRINTFUNCTION("Contains unrecognizable character\n");
           return -1;
        }
        score += pair2score[tmp];
    }
    return score;
}
double getES(vector<int>& indexes){
    double ES_max = -100000;
    double ES_current = 0;
    //sort(indexes.begin(), indexes.end());
    int N = score_matrix.size();
    double Nr = 0;
    int Nh = (int)predefined.size();
    for(int i=0; i < (int)indexes.size(); i++){
        if(score_matrix[indexes[i]].second > 0){
           Nr+=score_matrix[indexes[i]].second;
        }else{
           Nr-=score_matrix[indexes[i]].second;
        }
    }
    ES_current = (score_matrix[indexes[0]].second/Nr)-(double)indexes[0]/(N-Nh);
    //Log<<score_matrix[indexes[0]].first<<": "<<score_matrix[indexes[0]].second<<" ES"<<indexes[0]<<": "<<ES_current<<endl;
    if(ES_max < ES_current){
       ES_max = ES_current;
    }
    for(int i=1; i < (int)indexes.size(); i++){
        ES_current += (score_matrix[indexes[i]].second/Nr)-(double)(indexes[i]-indexes[i-1]-1)/(N-Nh);
        //Log<<score_matrix[indexes[i]].first<<": "<<score_matrix[indexes[i]].second<<" ES"<<indexes[i]<<": "<<ES_current<<endl;
        if(ES_max < ES_current){
           ES_max = ES_current;
        }
    }
    return ES_max;
}
double getES_curves(vector<int>& indexes){
    double ES_max = -100000;
    double ES_current = 0;
    //sort(indexes.begin(), indexes.end());
    int N = score_matrix.size();
    double Nr = 0;
    int Nh = (int)predefined.size();
    for(int i=0; i < (int)indexes.size(); i++){
        if(score_matrix[indexes[i]].second > 0){
           Nr+=score_matrix[indexes[i]].second;
        }else{
           Nr-=score_matrix[indexes[i]].second;
        }
    }
    line1 += "0\t";
    line2 += "0\t";
    ES_current = (score_matrix[indexes[0]].second/Nr)-(double)indexes[0]/(N-Nh);
    line1 += int2str(indexes[0])+"\t";
    line2 += double2str(ES_current-score_matrix[indexes[0]].second/Nr)+"\t";
    line1 += int2str(indexes[0]+1)+"\t";
    line2 += double2str(ES_current)+"\t";
    //Log<<score_matrix[indexes[0]].first<<": "<<score_matrix[indexes[0]].second<<" ES"<<indexes[0]<<": "<<ES_current<<endl;
    if(ES_max < ES_current){
       ES_max = ES_current;
    }
    for(int i=1; i < (int)indexes.size(); i++){
        ES_current += (score_matrix[indexes[i]].second/Nr)-(double)(indexes[i]-indexes[i-1]-1)/(N-Nh);
        //Log<<score_matrix[indexes[i]].first<<": "<<score_matrix[indexes[i]].second<<" ES"<<indexes[i]<<": "<<ES_current<<endl;
        if(ES_max < ES_current){
           ES_max = ES_current;
        }
        line1 += int2str(indexes[i])+"\t";
        line2 += double2str(ES_current-score_matrix[indexes[i]].second/Nr)+"\t";
        line1 += int2str(indexes[i]+1)+"\t";
        line2 += double2str(ES_current)+"\t";
    }
    line1 += int2str(score_matrix.size())+"\t";
    line2 += double2str(0)+"\t";
    //CURVE<<line1<<endl;
    //CURVE<<line2<<endl;
    return ES_max;
}
double getES_curves(){
    vector<int> indexes;
    for(int i=0; i < (int)score_matrix.size(); i++){
        if(predefined.count(score_matrix[i].first) != 0){
           indexes.push_back(i);
        }
    }
    //cerr<<"indexes: "<<indexes.size()<<endl;
    return getES_curves(indexes);
}
double getES(){
    vector<int> indexes;
    for(int i=0; i < (int)score_matrix.size(); i++){
        if(predefined.count(score_matrix[i].first) != 0){
           indexes.push_back(i);
        }
    }
    //cerr<<"indexes: "<<indexes.size()<<endl;
    return getES(indexes);
}
//vector<string> string_tokenize(const string& str, const string& delimiters = " \t\n\r", bool skip_empty = true);
inline vector<string> string_tokenize(const string& str, const string& delimiters, bool skip_empty) {
    // Skip delimiters at beginning.
        string::size_type lastPos = skip_empty ? str.find_first_not_of(delimiters, 0) : 0;
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
        vector<string> result;
        result.clear();

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
                //__ASSERT(pos > lastPos || !skip_empty, "internal error, pos <= lastPos.\n");

                //if (pos == lastPos) result.push_back("");
                result.push_back(str.substr(lastPos, pos - lastPos));

                if (pos == string::npos) break;
                if (pos == str.size() - 1) {
                        if (!skip_empty) result.push_back("");
                        break;
                }
        // Skip delimiters.  Note the "not_of"
                lastPos = skip_empty ? str.find_first_not_of(delimiters, pos) : pos + 1;
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
    return result;
}
string int2str(int value){
    char tmp[20];
    sprintf(tmp, "%d", value);
    return tmp;
}
string double2str(double value){
    char tmp[20];
    sprintf(tmp, "%.6f", value);
    return tmp;
}
void ToUpperString(string &str){
    transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper);
}
#ifdef INDEP_PROGRAM
#else
extern "C" {
int asebC(char ** input1, char ** input2, char ** input3, char ** output, int * Permuationtimes, int * proteins_or_sites)
{
    if(proteins_or_sites[0] == 0){
       Rprintf("background sites: %s\npredefined sites: %s\nsites to test: %s\noutput: %s\nPermutation times: %d\n", input1[0],input2[0],input3[0],output[0],Permuationtimes[0]);
    }else{
       Rprintf("background sites: %s\npredefined sites: %s\nproteins to test: %s\noutput: %s\nPermutation times: %d\n", input1[0],input2[0],input3[0],output[0],Permuationtimes[0]);
    }
    print_curves = 1;
    string curveFile = string(output[0])+".curves";
    CURVE.open(curveFile.c_str());
    if(!CURVE){
       Rprintf("Can not open %s\n", curveFile.c_str());
       return -1;
    }
    if(proteins_or_sites[0] == 0){
       aseb_sites(input1[0],input2[0],input3[0],output[0],Permuationtimes[0]);
    }else{
       aseb_protein(input1[0],input2[0],input3[0],output[0],Permuationtimes[0]);
    }
    pair2score.clear();
    geneName2rank.clear();
    predefined.clear();
    predefined_seq.clear();
    id2seq.clear();
    total_poteins.clear();
    p_times=1000;
    print_curves = 0;
    line1="";
    line2="";
    return 0;
}
}
extern "C" {
static R_NativePrimitiveArgType asebC_t[] = {STRSXP, STRSXP, STRSXP, STRSXP, INTSXP, INTSXP};
R_CMethodDef cMethods[] = {
  {".asebC", (DL_FUNC) &asebC, 6, asebC_t},
  {NULL, NULL, 0}
};

void R_init_ASEB(DllInfo *info){
   R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
}
#endif
