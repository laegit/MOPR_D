#include <iostream>
#include <map> 
#include <string>
#include <climits>
#include <stack>   
#include <fstream>      // std::ifstream 
#include <string>
#include <sstream>


using namespace std;

std::ofstream scriptsh;

template <typename T>
    std::string to_string(T value)
    {
      //create an output string stream
      std::ostringstream os ;

      //throw the value into the string stream
      os << value ;

      //convert the string stream into a string and return
      return os.str() ;
    }


void printRL(int n, int k, int i){
//int n = 30;
	//		int k = 3;
	//		int i = 1;
			for (int c = 1; c <= 35; ++c)
			{

				std::ofstream outt; 
				string nome = "st_"+to_string(n)+"_"+to_string(k)+"_"+to_string(i)+"_rl";
				scriptsh<<"sbatch "+nome<<endl;
				outt.open(nome, std::ofstream::out);
			
			outt<<"\#!/bin/bash"<<endl;
			outt<<"\#SBATCH --time=0-10:0   \# Especifica o tempo máximo de execução do job, dado no padrão dias-horas:minutos"<<endl; 
			
			
			outt<<"c="<<c<<endl;
			outt<<"n="<<n<<endl;
			outt<<"k="<<k<<endl;
			outt<<"i="<<i<<endl;
			//outt<<"mkdir Testes/VFSP/VFR_\"$j\"_\"$m\"/"<<endl;

			
			outt<<"ran=$(shuf -i 2000-9999999 -n 1);"<<endl;
			outt<<"./main_\"$k\"obj KCall/KC\"$n\"-\"$k\"fl-\"$i\"rl.dat $ran > Testes/KC\"$n\"_\"$k\"fl_\"$i\"rl/KC\"$n\"_\"$k\"fl_\"$i\"rl_\"$c\".out"<<endl;
			       
		
			}


}


void printuni(int n, int k, int i){
//int n = 30;
	//		int k = 3;
	//		int i = 1;
			for (int c = 1; c <= 35; ++c)
			{

				std::ofstream outt; 
				string nome = "st_"+to_string(n)+"_"+to_string(k)+"_"+to_string(i)+"_uni";
				scriptsh<<"sbatch "+nome<<endl;
				outt.open(nome, std::ofstream::out);
			
			outt<<"\#!/bin/bash"<<endl;
			outt<<"\#SBATCH --time=0-10:0   \# Especifica o tempo máximo de execução do job, dado no padrão dias-horas:minutos"<<endl; 
			
			
			outt<<"c="<<c<<endl;
			outt<<"n="<<n<<endl;
			outt<<"k="<<k<<endl;
			outt<<"i="<<i<<endl;
			//outt<<"mkdir Testes/VFSP/VFR_\"$j\"_\"$m\"/"<<endl;

			
			outt<<"ran=$(shuf -i 2000-9999999 -n 1);"<<endl;
			outt<<"./main_\"$k\"obj KCall/KC\"$n\"-\"$k\"fl-\"$i\"uni.dat $ran > Testes/KC\"$n\"_\"$k\"fl_\"$i\"uni/KC\"$n\"_\"$k\"fl_\"$i\"uni_\"$c\".out"<<endl;
			       
		
			}


}



int main( int argc,  char *argv[]){

			
			
		scriptsh.open("chamaTodosNVOO.sh", std::ofstream::out);
		printRL(30, 3, 1);
		printRL(30, 3, 2);
		printRL(30, 3, 3);
		printRL(40, 3, 1);
		printRL(40, 3, 2);
		printRL(40, 3, 3);
		printRL(50, 3, 1);
		printRL(50, 4, 1);
		printRL(50, 4, 2);
		printRL(50, 4, 3);
		printuni(30,3,1);
		printuni(30,3,2);
		printuni(30,3,3);
		printuni(40,3,1);
		printuni(40,3,2);
		printuni(40,3,3);
		printuni(50,4,1);
		printuni(50,4,2);
		printuni(50,4,3);


			
	

}